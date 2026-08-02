"""
analyze_stemness_patient_level.py
Reviewer-requested confirmatory analysis (manuscript comment #5): the
cell-level ACADL/CPT1A vs. stemness correlations in
analyze_acadl_stemness_correlation.py (r ~ 0.01-0.03 for CPT1A) risk
pseudoreplication -- p-values computed across tens of thousands of
non-independent cells from the same ~10-40 patients. This script re-tests
the same relationship two more rigorous ways per cohort:

  1. Patient-level pseudobulk: median CPT1A/ACADL expression and median
     stemness score per patient (malignant cells only), then Pearson/
     Spearman correlation ACROSS PATIENTS (n = number of patients, not
     number of cells).
  2. Mixed-effects model: gene ~ stemness_score + (1 | patient_id), fit on
     the cell-level data but with a random patient intercept, giving a
     p-value that accounts for within-patient non-independence directly
     (does not throw away cell-level resolution the way pseudobulk does).

Reuses the same local-scratch-copy convention as
analyze_acadl_stemness_correlation.py for the (large) h5ad files.
"""

import os
import shutil
import sys
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

np.random.seed(1234)

import scanpy as sc
import statsmodels.formula.api as smf
from scipy import stats
import yaml

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(BASE_DIR, "config")
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)
sc.settings.verbosity = 0

SCRATCH_DIR = os.environ.get(
    "PDAC_SCRATCH_DIR",
    r"C:\Users\erica\AppData\Local\Temp\claude\C--Users-erica-Onedrive-Research\2d87184b-c0c2-4c3e-8b0a-10233dd2fcb7\scratchpad\singlecell_scratch",
)
os.makedirs(SCRATCH_DIR, exist_ok=True)

FAO_GENES_OF_INTEREST = ["ACADL", "CPT1A"]
MIN_CELLS_PER_PATIENT = 30


def get_local_copy(processed_file_abs):
    basename = os.path.basename(processed_file_abs)
    local_path = os.path.join(SCRATCH_DIR, basename)
    if not os.path.exists(local_path):
        size_gb = os.path.getsize(processed_file_abs) / 1e9
        print(f"    Copying {basename} to local scratch ({size_gb:.2f} GB)...")
        shutil.copy2(processed_file_abs, local_path)
    else:
        print(f"    Local scratch copy already present: {local_path}")
    return local_path


def load_configs():
    with open(os.path.join(CONFIG_DIR, "singlecell_cohorts.yml")) as f:
        cohort_cfg = yaml.safe_load(f)
    with open(os.path.join(CONFIG_DIR, "gene_sets.yml")) as f:
        gene_sets = yaml.safe_load(f)
    return cohort_cfg, gene_sets


def analyze_cohort(cohort, stemness_genes):
    name = cohort["name"]
    processed_file = os.path.join(BASE_DIR, cohort["processed_file"])
    annotations_file = os.path.join(BASE_DIR, cohort["annotations_file"])

    if not os.path.exists(processed_file):
        print(f"  {name}: preprocessed file not found, skipping.")
        return []

    local_file = get_local_copy(processed_file)
    print(f"  Loading {name}...")
    adata = sc.read_h5ad(local_file)

    if os.path.exists(annotations_file):
        ann_df = pd.read_csv(annotations_file, sep="\t").set_index("cell_id")
        if "cell_type" not in adata.obs.columns:
            adata.obs["cell_type"] = "unknown"
        common = adata.obs.index.intersection(ann_df.index)
        adata.obs.loc[common, "cell_type"] = ann_df.loc[common, "cell_type"]

    mal = adata[adata.obs["cell_type"] == "malignant_epithelial"].copy()
    del adata

    stem_available = [g for g in stemness_genes if g in mal.var_names]
    if len(stem_available) == 0:
        print(f"  {name}: no stemness genes found, skipping.")
        return []

    np.random.seed(1234)
    sc.tl.score_genes(mal, stem_available, score_name="stemness_score", random_state=1234)

    counts = mal.obs["patient_id"].value_counts()
    keep_patients = counts[counts >= MIN_CELLS_PER_PATIENT].index
    keep_mask = mal.obs["patient_id"].isin(keep_patients).values
    print(f"  {name}: {keep_mask.sum()}/{mal.n_obs} cells retained "
          f"({len(keep_patients)} patients with >= {MIN_CELLS_PER_PATIENT} cells)")

    rows = []
    for gene in FAO_GENES_OF_INTEREST:
        if gene not in mal.var_names:
            print(f"    {gene}: not found in {name}, skipping.")
            continue
        expr = mal[:, gene].X
        expr = np.asarray(expr.todense()).flatten() if hasattr(expr, "todense") else np.asarray(expr).flatten()

        cell_df = pd.DataFrame({
            "expr": expr,
            "stemness_score": mal.obs["stemness_score"].values,
            "patient_id": mal.obs["patient_id"].values,
        })
        cell_df_f = cell_df[cell_df["patient_id"].isin(keep_patients)]

        if len(keep_patients) < 3:
            print(f"    {gene}: fewer than 3 patients with enough cells, skipping patient-level tests.")
            continue

        # 1. Patient-level pseudobulk correlation. ACADL/CPT1A are sparsely
        # detected at the single-cell level (unlike multi-gene signature
        # scores), so the per-patient MEDIAN of a single gene is 0 (zero
        # variance) in most cohorts; use the mean instead.
        pb = cell_df_f.groupby("patient_id", observed=True).agg(
            expr_mean=("expr", "mean"),
            stemness_mean=("stemness_score", "mean"),
            n_cells=("expr", "size"),
        )
        pb = pb[pb["n_cells"] > 0]
        pb_r, pb_p = stats.pearsonr(pb["expr_mean"], pb["stemness_mean"])
        pb_rho, pb_rho_p = stats.spearmanr(pb["expr_mean"], pb["stemness_mean"])

        # 2. Mixed-effects model with patient random intercept
        try:
            md = smf.mixedlm("expr ~ stemness_score", cell_df_f, groups=cell_df_f["patient_id"])
            mdf = md.fit(reml=True)
            me_coef = float(mdf.params["stemness_score"])
            me_p = float(mdf.pvalues["stemness_score"])
            ci = mdf.conf_int().loc["stemness_score"]
            me_ci_lo, me_ci_hi = float(ci[0]), float(ci[1])
        except Exception as e:
            print(f"    {gene}: mixed model failed ({e}), leaving blank.")
            me_coef = me_p = me_ci_lo = me_ci_hi = np.nan

        print(f"    {gene}: patient-pseudobulk Pearson r={pb_r:.4f} (p={pb_p:.4g}, n={len(pb)} patients); "
              f"mixed-model coef={me_coef:.5f} (p={me_p:.4g})")

        rows.append({
            "cohort": name, "gene": gene, "n_patients": len(pb), "n_cells_used": int(pb["n_cells"].sum()),
            "patient_pseudobulk_pearson_r": round(pb_r, 4), "patient_pseudobulk_pearson_p": pb_p,
            "patient_pseudobulk_spearman_rho": round(pb_rho, 4), "patient_pseudobulk_spearman_p": pb_rho_p,
            "mixedlm_stemness_coef": round(me_coef, 6) if not np.isnan(me_coef) else np.nan,
            "mixedlm_stemness_ci_lower": round(me_ci_lo, 6) if not np.isnan(me_ci_lo) else np.nan,
            "mixedlm_stemness_ci_upper": round(me_ci_hi, 6) if not np.isnan(me_ci_hi) else np.nan,
            "mixedlm_stemness_pvalue": me_p,
        })

    del mal
    return rows


def main():
    print("=== ACADL/CPT1A vs Stemness: Patient-Level Pseudobulk + Mixed-Effects Model ===\n")
    cohort_cfg, gene_sets = load_configs()
    stemness_genes = gene_sets["pdac_stemness"]

    all_rows = []
    for cohort in cohort_cfg["cohorts"]:
        print(f"--- {cohort['name']} ---")
        try:
            all_rows.extend(analyze_cohort(cohort, stemness_genes))
        except Exception as e:
            import traceback
            print(f"  ERROR in {cohort['name']}: {e}")
            traceback.print_exc()
        print()

    if not all_rows:
        print("No results produced.")
        return

    out = pd.DataFrame(all_rows)
    out_path = os.path.join(TABLES_DIR, "acadl_cpt1a_stemness_patient_level.tsv")
    out.to_csv(out_path, sep="\t", index=False)
    print(f"Results saved: {out_path}")
    print(out.to_string(index=False))
    print("\n=== Complete ===")


if __name__ == "__main__":
    main()
