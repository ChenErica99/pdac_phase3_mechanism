"""
analyze_acadl_malignant_pseudobulk.py
Reviewer-requested confirmatory analysis (manuscript comment #3): test
whether the bulk/CPTAC ACADL reduction in "aggressive" (hypoxia-high/
acinar-low) tumors is tumor-intrinsic, by testing ACADL (a) across cell
types, and (b) in patient-level malignant-cell-ONLY pseudobulk, using the
three real single-cell/single-nucleus cohorts already used elsewhere in
this pipeline (Section 2.2-2.3).

This cannot use the CPTAC cohort directly (no matched single-cell data for
those patients); instead it asks the same acinar-content-confound question
within the single-cell cohorts, where a patient's own hypoxia/acinar group
label and cell-type composition are both directly observable:

  (a) Across cell types: is ACADL simply an acinar-lineage marker (much
      higher in acinar_normal cells than malignant cells)? If so, a
      composite bulk sample with less acinar content would show lower
      ACADL for a purely compositional reason.
  (b) Malignant-only pseudobulk: restricting entirely to malignant cells
      (excluding acinar_normal, ductal_normal, stroma, immune), does
      ACADL pseudobulk still differ between patients whose OWN malignant
      cells are hypoxia-high/acinar-low ("aggressive") vs hypoxia-low/
      acinar-high ("reference")? A positive answer here is direct
      evidence the ACADL effect is not purely acinar-content-driven.

Patient group labels are derived from each patient's own malignant-cell
pseudobulk hypoxia_score/acinar_identity_score (median split within
cohort), reusing the same MIN_CELLS_PER_PATIENT/MIN_PATIENTS_PER_ARM
convention as analyze_lipid_cell_of_origin.py.
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
from scipy import stats
from statsmodels.stats.multitest import multipletests
import yaml

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(BASE_DIR, "config")
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
SC_DIR = os.path.join(BASE_DIR, "data", "processed", "singlecell")
os.makedirs(TABLES_DIR, exist_ok=True)
sc.settings.verbosity = 0

MIN_CELLS_PER_PATIENT = 5
MIN_PATIENTS_PER_ARM = 3
GENE = "ACADL"

SCRATCH_DIR = os.environ.get(
    "PDAC_SCRATCH_DIR",
    r"C:\Users\erica\AppData\Local\Temp\claude\C--Users-erica-Onedrive-Research\2d87184b-c0c2-4c3e-8b0a-10233dd2fcb7\scratchpad\singlecell_scratch",
)
os.makedirs(SCRATCH_DIR, exist_ok=True)


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
        return yaml.safe_load(f)["cohorts"]


def analyze_cohort(cohort):
    name = cohort["name"]
    processed_file = os.path.join(BASE_DIR, cohort["processed_file"])
    annotations_file = os.path.join(BASE_DIR, cohort["annotations_file"])
    scores_file = os.path.join(BASE_DIR, cohort["scores_file"])

    if not (os.path.exists(processed_file) and os.path.exists(scores_file)):
        print(f"  {name}: missing processed/scores file, skipping.")
        return None, None

    local_file = get_local_copy(processed_file)
    print(f"  Loading {name}...")
    adata = sc.read_h5ad(local_file)

    if GENE not in adata.var_names:
        print(f"  {name}: {GENE} not in var_names, skipping.")
        del adata
        return None, None

    scores = pd.read_csv(scores_file, sep="\t").set_index("cell_id")
    if os.path.exists(annotations_file):
        ann = pd.read_csv(annotations_file, sep="\t").set_index("cell_id")
        if "cell_type" not in scores.columns and "cell_type" in ann.columns:
            scores["cell_type"] = ann["cell_type"].reindex(scores.index)

    expr = adata[:, GENE].X
    expr = np.asarray(expr.todense()).flatten() if hasattr(expr, "todense") else np.asarray(expr).flatten()
    cell_df = pd.DataFrame({"ACADL": expr}, index=adata.obs.index)
    cell_df = cell_df.join(scores[["patient_id", "cell_type", "hypoxia_score", "acinar_identity_score"]], how="inner")
    del adata

    # ---- (a) ACADL across cell types (all patients pooled) ----
    ct_summary = (
        cell_df.groupby("cell_type", observed=True)["ACADL"]
        .agg(["mean", "median", "count"])
        .reset_index()
        .rename(columns={"mean": "acadl_mean", "median": "acadl_median", "count": "n_cells"})
    )
    ct_summary.insert(0, "cohort", name)

    # ---- (b) malignant-only patient pseudobulk, grouped by patient's OWN
    # malignant-cell hypoxia/acinar status ----
    mal = cell_df[cell_df["cell_type"] == "malignant_epithelial"]
    patient_n = mal.groupby("patient_id").size()
    valid_patients = patient_n[patient_n >= MIN_CELLS_PER_PATIENT].index
    mal = mal[mal["patient_id"].isin(valid_patients)]

    # ACADL is sparsely detected at the single-cell level (unlike the
    # multi-gene signature scores used elsewhere in this pipeline), so a
    # per-patient MEDIAN of a single gene is 0 in most patients/cohorts
    # and uninformative; use the mean instead, which is sensitive to the
    # (real, non-zero) fraction of expressing cells.
    pseudobulk = mal.groupby("patient_id", observed=True).agg(
        acadl_pseudobulk=("ACADL", "mean"),
        hypoxia_score=("hypoxia_score", "median"),
        acinar_score=("acinar_identity_score", "median"),
        n_malignant_cells=("ACADL", "size"),
    )

    if len(pseudobulk) < 2 * MIN_PATIENTS_PER_ARM:
        print(f"  {name}: only {len(pseudobulk)} patients with >={MIN_CELLS_PER_PATIENT} malignant cells, "
              f"insufficient for group comparison.")
        return ct_summary, None

    hyp_med = pseudobulk["hypoxia_score"].median()
    acin_med = pseudobulk["acinar_score"].median()
    hyp_high = pseudobulk["hypoxia_score"] >= hyp_med
    acin_low = pseudobulk["acinar_score"] < acin_med
    pseudobulk["group"] = "other"
    pseudobulk.loc[hyp_high & acin_low, "group"] = "aggressive"
    pseudobulk.loc[~hyp_high & ~acin_low, "group"] = "reference"

    agg = pseudobulk.loc[pseudobulk["group"] == "aggressive", "acadl_pseudobulk"]
    ref = pseudobulk.loc[pseudobulk["group"] == "reference", "acadl_pseudobulk"]
    print(f"  {name}: malignant-only pseudobulk — aggressive n={len(agg)} patients, reference n={len(ref)} patients")

    result_row = None
    if len(agg) >= MIN_PATIENTS_PER_ARM and len(ref) >= MIN_PATIENTS_PER_ARM:
        stat, p = stats.ranksums(agg, ref)
        direction = "down" if agg.median() < ref.median() else "up"
        print(f"    ACADL malignant-only pseudobulk: aggressive median={agg.median():.4f}, "
              f"reference median={ref.median():.4f}, direction={direction}, p={p:.4g}")
        result_row = {
            "cohort": name, "n_aggressive_patients": len(agg), "n_reference_patients": len(ref),
            "acadl_median_aggressive": round(float(agg.median()), 4),
            "acadl_median_reference": round(float(ref.median()), 4),
            "direction": direction, "p_value": round(float(p), 6),
        }
    else:
        print(f"    {name}: fewer than {MIN_PATIENTS_PER_ARM} patients in one arm, skipping test "
              f"(aggressive={len(agg)}, reference={len(ref)}).")

    return ct_summary, result_row


def main():
    print("=== ACADL Across Cell Types + Malignant-Only Patient Pseudobulk (single-cell cohorts) ===\n")
    cohorts = load_configs()

    ct_rows, pb_rows = [], []
    for cohort in cohorts:
        print(f"--- {cohort['name']} ---")
        try:
            ct_summary, pb_row = analyze_cohort(cohort)
            if ct_summary is not None:
                ct_rows.append(ct_summary)
            if pb_row is not None:
                pb_rows.append(pb_row)
        except Exception as e:
            import traceback
            print(f"  ERROR in {cohort['name']}: {e}")
            traceback.print_exc()
        print()

    if ct_rows:
        ct_df = pd.concat(ct_rows, ignore_index=True)
        ct_out = os.path.join(TABLES_DIR, "acadl_by_celltype.tsv")
        ct_df.to_csv(ct_out, sep="\t", index=False)
        print(f"ACADL-by-cell-type table saved: {ct_out}")
        print(ct_df.to_string(index=False))

    if pb_rows:
        pb_df = pd.DataFrame(pb_rows)
        if len(pb_df) > 0:
            _, fdr, _, _ = multipletests(pb_df["p_value"].values, method="fdr_bh")
            pb_df["fdr"] = fdr.round(4)
        pb_out = os.path.join(TABLES_DIR, "acadl_malignant_only_pseudobulk.tsv")
        pb_df.to_csv(pb_out, sep="\t", index=False)
        print(f"\nMalignant-only pseudobulk table saved: {pb_out}")
        print(pb_df.to_string(index=False))
    else:
        print("\nNo cohort had sufficient patients per arm for the malignant-only pseudobulk test.")

    print("\n=== Complete ===")


if __name__ == "__main__":
    main()
