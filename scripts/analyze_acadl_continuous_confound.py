"""
analyze_acadl_continuous_confound.py
Reviewer-requested confirmatory analysis (manuscript comment #3, most
important biological confounding issue): the aggressive subgroup is
partly DEFINED by low acinar-identity score, so lower bulk ACADL in that
subgroup could partly reflect reduced normal-acinar tissue content /
lineage composition rather than a tumor-intrinsic change.

Fits, per proteomic source (umich, bcm) and separately for ACADL protein
and ACADL mRNA:

    gene ~ hypoxia_z + acinar_z + hypoxia_z:acinar_z + purity + grade + stage

using ALL CPTAC tumor samples with a valid continuous hypoxia/acinar score
(not just the discrete aggressive/reference median-split groups), so the
independent contribution of the acinar axis can be estimated directly
rather than folded into a single group indicator. The median-split
aggressive-vs-reference comparison (analyze_cptac_protein.py) is retained
as a sensitivity analysis, per the reviewer's request, not replaced.

Reuses load_cptac_dataset/load_proteomics/zscore_genes from
analyze_cptac_protein.py rather than reimplementing CPTAC loading.
"""

import os
import sys
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

np.random.seed(1234)

import statsmodels.api as sm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_cptac_protein import (
    load_cptac_dataset, load_proteomics, zscore_genes,
    HYPOXIA_GENES, ACINAR_GENES,
)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

GENE = "ACADL"


def get_gene_series(df, gene):
    col = df[gene]
    if isinstance(col, pd.DataFrame):
        col = col.mean(axis=1)
    return col


def continuous_scores(trans):
    hyp_avail = [g for g in HYPOXIA_GENES if g in trans.columns]
    acin_avail = [g for g in ACINAR_GENES if g in trans.columns]
    hyp_z = zscore_genes(trans[hyp_avail]).mean(axis=1)
    acin_z = zscore_genes(trans[acin_avail]).mean(axis=1)
    return pd.DataFrame({"hypoxia_z": hyp_z, "acinar_z": acin_z})


def fit_model(y, scores, purity_df, clin, label):
    common = y.index.intersection(scores.index)
    if purity_df is not None:
        common = common.intersection(purity_df.index)
    if clin is not None:
        common = common.intersection(clin.index)

    df = pd.DataFrame(index=common)
    df["y"] = y.loc[common].astype(float)
    df["hypoxia_z"] = scores.loc[common, "hypoxia_z"]
    df["acinar_z"] = scores.loc[common, "acinar_z"]
    df["interaction"] = df["hypoxia_z"] * df["acinar_z"]
    df["purity"] = purity_df.loc[common, "TumorPurity"] if purity_df is not None else np.nan
    df["grade"] = clin.loc[common, "grade_numeric"] if clin is not None else np.nan
    df["stage"] = clin.loc[common, "stage_numeric"] if clin is not None else np.nan

    df = df.dropna(subset=["y"])
    for c in ["purity", "grade", "stage"]:
        if df[c].notna().any():
            df[c] = df[c].fillna(df[c].median())
    df = df.dropna()

    n = len(df)
    print(f"  [{label}] n={n} samples with complete data")
    if n < 15:
        print(f"  [{label}] too few samples, skipping.")
        return None

    X = sm.add_constant(df[["hypoxia_z", "acinar_z", "interaction", "purity", "grade", "stage"]])
    model = sm.OLS(df["y"], X).fit()

    row = {"label": label, "n": n}
    for term in ["hypoxia_z", "acinar_z", "interaction", "purity", "grade", "stage"]:
        row[f"coef_{term}"] = round(float(model.params[term]), 5)
        row[f"pval_{term}"] = round(float(model.pvalues[term]), 6)
    row["r_squared"] = round(float(model.rsquared), 4)
    print(f"    hypoxia_z: coef={row['coef_hypoxia_z']:.4f} p={row['pval_hypoxia_z']:.4g}  "
          f"acinar_z: coef={row['coef_acinar_z']:.4f} p={row['pval_acinar_z']:.4g}  "
          f"interaction: coef={row['coef_interaction']:.4f} p={row['pval_interaction']:.4g}")
    return row, model


def main():
    print("=== ACADL Continuous Confound Model: hypoxia + acinar + interaction + purity + grade + stage ===\n")
    pdac, trans, purity_df, clin = load_cptac_dataset()

    if GENE not in trans.columns:
        print(f"{GENE} not in washu transcriptomics, aborting.")
        return

    scores = continuous_scores(trans)
    rna = get_gene_series(trans, GENE)

    rows = []
    for source in ["umich", "bcm"]:
        print(f"\n--- {source} protein ---")
        prot = load_proteomics(pdac, source)
        if GENE in prot.columns:
            result = fit_model(prot[GENE], scores, purity_df, clin, f"{source}_protein")
            if result is not None:
                rows.append(result[0])
        else:
            print(f"  {GENE} not in {source} proteomics.")

    print("\n--- washu RNA ---")
    result = fit_model(rna, scores, purity_df, clin, "washu_rna")
    if result is not None:
        rows.append(result[0])

    if not rows:
        print("\nNo models fit. Exiting.")
        return

    out = pd.DataFrame(rows)
    out_path = os.path.join(TABLES_DIR, "acadl_continuous_confound_model.tsv")
    out.to_csv(out_path, sep="\t", index=False)
    print(f"\nResults saved: {out_path}")
    print(out.to_string(index=False))

    print("\nInterpretation guide: a significant, negative-going hypoxia_z coefficient with a "
          "non-significant acinar_z coefficient (after mutual adjustment) would support hypoxia, "
          "not acinar-lineage content, as the operative axis for ACADL suppression. A significant "
          "positive acinar_z coefficient (higher acinar score -> higher ACADL) with attenuated "
          "hypoxia_z would instead support the acinar-lineage-content confound the reviewer raised.")

    print("\n=== Complete ===")


if __name__ == "__main__":
    main()
