#!/usr/bin/env python3
"""Cross-DATASET-validated gene panel on each stratum's meta signature.

Honest validation: Leave-One-Dataset-Out CV (a whole study is held out each fold),
NOT leave-one-sample-out. RNA-seq and microarray are kept separate (different scales).
Per-dataset z-scoring removes between-study scale/batch before combining.

Inputs (produced by the pipeline):
  results/meta/<stratum>.meta.tsv        meta signature (gene, meta_p/meta_padj, ...)
  data/expr/<dataset>.expr.csv.gz        genes x samples (from fetch_expr.R)
  data/expr/<dataset>.meta.csv           samples x covariates (from fetch_expr.R)
  config/samplesheet.csv                 dataset -> condition/tissue/technique/group labels

Output:
  results/ml/panel_report.tsv            one row per (stratum, technique)
8 GB note: features are restricted to the meta-signature genes, so matrices stay small.
"""
import argparse, glob, os, sys, warnings
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")

def norm_id(cond, tissue):
    return f"{cond}__{tissue}".replace(" ", "_").replace("/", "-")

def signature_genes(meta_path, padj, max_genes):
    m = pd.read_csv(meta_path, sep="\t")
    if "gene" not in m.columns or len(m) == 0:
        return []
    col = "meta_padj" if "meta_padj" in m.columns else ("min_padj" if "min_padj" in m.columns else "meta_p")
    if col in m.columns:
        m = m[pd.to_numeric(m[col], errors="coerce").notna()]
        m = m.sort_values(col)
        sig = m[pd.to_numeric(m[col], errors="coerce") < padj]
        genes = (sig if len(sig) > 0 else m).head(max_genes)["gene"].astype(str).tolist()
    else:
        genes = m["gene"].astype(str).head(max_genes).tolist()
    return list(dict.fromkeys(genes))

def load_dataset(expr_dir, ds, genes):
    ef = os.path.join(expr_dir, f"{ds}.expr.csv.gz")
    mf = os.path.join(expr_dir, f"{ds}.meta.csv")
    if not (os.path.exists(ef) and os.path.exists(mf)):
        return None, None
    expr = pd.read_csv(ef, index_col=0)              # genes x samples
    expr = expr[~expr.index.duplicated(keep="first")]
    keep = [g for g in genes if g in expr.index]
    if not keep:
        return None, None
    X = expr.loc[keep].T                              # samples x genes
    meta = pd.read_csv(mf, index_col=0)               # samples x covariates
    return X, meta

def zscore_per_dataset(X):
    mu = X.mean(axis=0); sd = X.std(axis=0).replace(0, np.nan)
    Z = (X - mu) / sd
    return Z.fillna(0.0)

def run_stratum(stratum, ss, expr_dir, genes, min_datasets):
    rows = []
    sub = ss[ss["__sid"] == stratum]
    for tech, g in sub.groupby("technique"):
        Xs, ys, grps = [], [], []
        used = 0
        for _, r in g.iterrows():
            ds = str(r["dataset_id"])
            X, meta = load_dataset(expr_dir, ds, genes)
            if X is None:
                continue
            gc, ca, co = r.get("group_column"), str(r.get("case_label")), str(r.get("control_label"))
            if pd.isna(gc) or gc not in meta.columns:
                continue
            lab = meta[gc].astype(str)
            y = pd.Series(index=X.index, dtype="float")
            common = X.index.intersection(meta.index)
            X = X.loc[common]; lab = lab.loc[common]
            y = lab.map({ca: 1, co: 0})
            mask = y.notna()
            if mask.sum() < 4 or y[mask].nunique() < 2:
                continue
            Xz = zscore_per_dataset(X[mask])
            Xs.append(Xz); ys.append(y[mask].astype(int)); grps += [ds] * int(mask.sum()); used += 1
        if used < max(2, min_datasets // 2):   # need >=2 datasets of this technique to cross-validate
            rows.append(dict(stratum=stratum, technique=tech, n_datasets=used, status="skipped (<2 usable datasets)"))
            continue
        Xall = pd.concat(Xs, axis=0)
        Xall = Xall.reindex(columns=sorted(set().union(*[set(x.columns) for x in Xs])))
        yall = pd.concat(ys).values
        groups = np.array(grps)
        logo = LeaveOneGroupOut()
        clf = Pipeline([("imp", SimpleImputer(strategy="mean")),
                        ("sc", StandardScaler()),
                        ("lr", LogisticRegression(penalty="elasticnet", solver="saga",
                                                  l1_ratio=0.5, C=1.0, max_iter=5000))])
        aucs = []
        for tr, te in logo.split(Xall.values, yall, groups):
            if len(np.unique(yall[te])) < 2 or len(np.unique(yall[tr])) < 2:
                continue
            clf.fit(Xall.values[tr], yall[tr])
            p = clf.predict_proba(Xall.values[te])[:, 1]
            aucs.append(roc_auc_score(yall[te], p))
        # final panel on all data
        clf.fit(Xall.values, yall)
        coef = clf.named_steps["lr"].coef_.ravel()
        panel = pd.Series(np.abs(coef), index=Xall.columns)
        nonzero = panel[panel > 0].sort_values(ascending=False)
        rows.append(dict(
            stratum=stratum, technique=tech, n_datasets=used, n_samples=len(yall),
            n_signature_genes=Xall.shape[1],
            mean_auroc=round(float(np.mean(aucs)), 3) if aucs else np.nan,
            sd_auroc=round(float(np.std(aucs)), 3) if aucs else np.nan,
            n_folds_scored=len(aucs), panel_size=int((panel > 0).sum()),
            top_panel_genes=";".join(nonzero.head(25).index.tolist()),
            status="ok"))
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta_dir", default="results/meta")
    ap.add_argument("--expr_dir", default="data/expr")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--cv", default="leave-one-dataset-out")
    ap.add_argument("--padj", type=float, default=0.05)
    ap.add_argument("--max_genes", type=int, default=300)
    ap.add_argument("--min_datasets", type=int, default=4)
    ap.add_argument("--out", default="results/ml/panel_report.tsv")
    a = ap.parse_args()

    ss = pd.read_csv(a.samplesheet)
    if "notes" in ss.columns:
        ss = ss[~ss["notes"].astype(str).str.contains("DROP", case=False, na=False)]
    ss["__sid"] = [norm_id(c, t) for c, t in zip(ss["condition"], ss["tissue_primary"])]

    rows = []
    for mp in sorted(glob.glob(os.path.join(a.meta_dir, "*.meta.tsv"))):
        stratum = os.path.basename(mp)[:-len(".meta.tsv")]
        genes = signature_genes(mp, a.padj, a.max_genes)
        if not genes:
            rows.append(dict(stratum=stratum, status="no signature genes")); continue
        rows += run_stratum(stratum, ss, a.expr_dir, genes, a.min_datasets)

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"\nWrote {a.out}")

if __name__ == "__main__":
    main()
