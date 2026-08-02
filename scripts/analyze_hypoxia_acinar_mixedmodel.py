"""
analyze_hypoxia_acinar_mixedmodel.py
Reviewer-requested confirmatory analysis (manuscript comment #2): a formal
mixed-effects model for the hypoxia-acinar single-cell co-occurrence test
(Section 2.2), replacing/augmenting the fixed-effect average of per-patient
Pearson correlations with a proper random-intercept mixed model:

    acinar_identity_score ~ hypoxia_score + (1 | patient_id)

fit separately per cohort on malignant cells only, using the already-
computed per-cell signature scores (*_cell_scores.tsv) -- no need to
reload the h5ad objects.
"""

import os
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

import statsmodels.formula.api as smf
from scipy import stats

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC_DIR = os.path.join(BASE_DIR, "data", "processed", "singlecell")
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

COHORTS = ["GSE154778", "GSE202051", "Peng_et_al"]
MIN_CELLS_PER_PATIENT = 30


def analyze_cohort(name):
    scores_file = os.path.join(SC_DIR, f"{name}_cell_scores.tsv")
    if not os.path.exists(scores_file):
        print(f"  {name}: scores file not found, skipping.")
        return None

    df = pd.read_csv(scores_file, sep="\t")
    mal = df[df["cell_type"] == "malignant_epithelial"].copy()
    counts = mal["patient_id"].value_counts()
    keep = counts[counts >= MIN_CELLS_PER_PATIENT].index
    mal = mal[mal["patient_id"].isin(keep)]
    n_patients = mal["patient_id"].nunique()
    print(f"  {name}: {len(mal)} malignant cells across {n_patients} patients "
          f"(>= {MIN_CELLS_PER_PATIENT} cells/patient)")

    if n_patients < 3:
        print(f"  {name}: too few patients for a mixed model, skipping.")
        return None

    # Pooled (naive) correlation, for reference against the mixed model
    pooled_r, pooled_p = stats.pearsonr(mal["hypoxia_score"], mal["acinar_identity_score"])

    md = smf.mixedlm(
        "acinar_identity_score ~ hypoxia_score",
        mal,
        groups=mal["patient_id"],
    )
    mdf = md.fit(reml=True)
    coef = mdf.params["hypoxia_score"]
    pval = mdf.pvalues["hypoxia_score"]
    ci_lo, ci_hi = mdf.conf_int().loc["hypoxia_score"]
    re_var = float(mdf.cov_re.iloc[0, 0])
    resid_var = float(mdf.scale)
    icc = re_var / (re_var + resid_var) if (re_var + resid_var) > 0 else np.nan

    print(f"    Mixed model (random patient intercept): hypoxia_score coef={coef:.4f} "
          f"[{ci_lo:.4f}, {ci_hi:.4f}], p={pval:.4g}; ICC(patient)={icc:.3f}")
    print(f"    Pooled cell-level Pearson r={pooled_r:.4f} (p={pooled_p:.4g}), for comparison")

    return {
        "cohort": name,
        "n_cells": len(mal),
        "n_patients": n_patients,
        "pooled_pearson_r": round(pooled_r, 4),
        "pooled_pearson_p": pooled_p,
        "mixedlm_hypoxia_coef": round(float(coef), 5),
        "mixedlm_hypoxia_ci_lower": round(float(ci_lo), 5),
        "mixedlm_hypoxia_ci_upper": round(float(ci_hi), 5),
        "mixedlm_hypoxia_pvalue": float(pval),
        "icc_patient": round(float(icc), 4),
    }


def main():
    print("=== Hypoxia-Acinar Mixed-Effects Model (patient random intercept, malignant cells) ===\n")
    rows = []
    for name in COHORTS:
        print(f"--- {name} ---")
        row = analyze_cohort(name)
        if row is not None:
            rows.append(row)
        print()

    if not rows:
        print("No cohorts produced results.")
        return

    out = pd.DataFrame(rows)
    out_path = os.path.join(TABLES_DIR, "hypoxia_acinar_mixedmodel.tsv")
    out.to_csv(out_path, sep="\t", index=False)
    print(f"Results saved: {out_path}")
    print(out.to_string(index=False))
    print("\n=== Complete ===")


if __name__ == "__main__":
    main()
