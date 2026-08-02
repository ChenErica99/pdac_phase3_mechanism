"""
analyze_acadl_stemness_correlation.py
Tests whether ACADL suppression in bulk/aggressive tumors (Sections 3-7) is
compatible with the literature claim that fatty acid oxidation is required
for pancreatic cancer stem cell (PaCSC) tumorigenicity (PMC11351511, Fatty
acid oxidation is critical for the tumorigenic potential and chemoresistance
of pancreatic cancer stem cells, J Transl Med 2024).

That paper's core finding is CPT1A-centered (not ACADL) and is measured in a
CD133+/CD44+/sphere-forming stem-cell-ENRICHED subpopulation via PDX-derived
cultures and circulating tumor cells, not bulk tumor tissue. This pipeline's
ACADL finding (Sections 3-7) is measured in bulk tumor proteomics/RNA-seq
(CPTAC) and averaged single-cell malignant compartments — dominated by the
non-stem majority. A real reconciliation, not just a caveat, requires testing
whether ACADL (and CPT1A, the paper's own gene) actually goes UP with
increasing stemness within the malignant-cell compartment, i.e. whether a
bulk-tissue-average finding and a stem-subpopulation finding can both be true
because they describe different cell fractions of the same tumor.

Method: score each malignant cell's stemness using the SAME marker panel the
literature paper uses to define its CD133+/CD44+ CSC population (PROM1=CD133,
CD44, NANOG, KLF4, POU5F1=OCT3/4, SOX2 — see config/gene_sets.yml
`pdac_stemness`), then correlate raw ACADL (and CPT1A) expression against
that score across all malignant cells in each of the 3 real scRNA-seq
cohorts. This is a cell-level correlation (like Section 3's hypoxia/acinar
co-occurrence test), not a patient-arm group comparison, so it is NOT subject
to the pseudoreplication issue that affected the Section 4/5 cell-of-origin
tests -- no patient-level aggregation or MIN_PATIENTS_PER_ARM threshold
applies here.
"""

import os
import shutil
import sys
import numpy as np
import pandas as pd
import random
import yaml
import warnings
warnings.filterwarnings("ignore")

np.random.seed(1234)
random.seed(1234)

import scanpy as sc
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(BASE_DIR, "config")
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
FIGURES_DIR = os.path.join(BASE_DIR, "results", "figures")
os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)
sc.settings.verbosity = 0

# Local (non-OneDrive) scratch dir — GSE202051's h5ad is ~5GB and reading it
# repeatedly from inside the OneDrive-synced project folder is catastrophically
# slow (see config/singlecell_cohorts.yml comment and Phase 3 memory). Copy
# each h5ad here before reading, delete when done, same convention as
# sensitivity_treatment_status_GSE202051.py used.
SCRATCH_DIR = os.environ.get(
    "PDAC_SCRATCH_DIR",
    r"C:\Users\erica\AppData\Local\Temp\claude\C--Users-erica-OneDrive-Research\753e4173-3eb7-4967-8750-b24a4cb9fc93\scratchpad\singlecell_scratch",
)
os.makedirs(SCRATCH_DIR, exist_ok=True)

FAO_GENES_OF_INTEREST = ["ACADL", "CPT1A"]  # ACADL = this pipeline's anchor; CPT1A = the literature paper's gene


def load_configs():
    with open(os.path.join(CONFIG_DIR, "singlecell_cohorts.yml")) as f:
        cohort_cfg = yaml.safe_load(f)
    with open(os.path.join(CONFIG_DIR, "gene_sets.yml")) as f:
        gene_sets = yaml.safe_load(f)
    return cohort_cfg, gene_sets


def get_local_copy(processed_file_abs):
    """Copy the h5ad to local scratch if not already there; return local path."""
    basename = os.path.basename(processed_file_abs)
    local_path = os.path.join(SCRATCH_DIR, basename)
    if not os.path.exists(local_path):
        size_gb = os.path.getsize(processed_file_abs) / 1e9
        print(f"    Copying {basename} to local scratch ({size_gb:.2f} GB)...")
        shutil.copy2(processed_file_abs, local_path)
    else:
        print(f"    Local scratch copy already present: {local_path}")
    return local_path


def analyze_cohort(cohort, stemness_genes):
    name = cohort["name"]
    processed_file = os.path.join(BASE_DIR, cohort["processed_file"])
    annotations_file = os.path.join(BASE_DIR, cohort["annotations_file"])

    if not os.path.exists(processed_file):
        print(f"  {name}: preprocessed file not found, skipping.")
        return None, None

    local_file = get_local_copy(processed_file)
    print(f"  Loading {name} from local scratch...")
    adata = sc.read_h5ad(local_file)

    if os.path.exists(annotations_file):
        ann_df = pd.read_csv(annotations_file, sep="\t").set_index("cell_id")
        if "cell_type" not in adata.obs.columns:
            adata.obs["cell_type"] = "unknown"
        common = adata.obs.index.intersection(ann_df.index)
        adata.obs.loc[common, "cell_type"] = ann_df.loc[common, "cell_type"]

    mal = adata[adata.obs["cell_type"] == "malignant_epithelial"].copy()
    n_patients = mal.obs["patient_id"].nunique() if "patient_id" in mal.obs.columns else np.nan
    print(f"  {name}: {mal.n_obs} malignant cells, {n_patients} patients")

    stem_available = [g for g in stemness_genes if g in mal.var_names]
    print(f"  {name}: stemness genes found {len(stem_available)}/{len(stemness_genes)}: {stem_available}")
    if len(stem_available) == 0:
        print(f"  {name}: no stemness marker genes found, skipping.")
        del adata, mal
        return None, None

    np.random.seed(1234)
    sc.tl.score_genes(mal, stem_available, score_name="stemness_score", random_state=1234)
    stemness = mal.obs["stemness_score"].values

    rows = []
    plot_data = {"stemness_score": stemness}
    for gene in FAO_GENES_OF_INTEREST:
        if gene not in mal.var_names:
            print(f"    {gene}: not found in {name}, skipping.")
            continue
        expr = mal[:, gene].X
        expr = np.asarray(expr.todense()).flatten() if hasattr(expr, "todense") else np.asarray(expr).flatten()
        plot_data[gene] = expr

        pearson_r, pearson_p = stats.pearsonr(expr, stemness)
        spearman_r, spearman_p = stats.spearmanr(expr, stemness)

        # Top-decile vs bottom-decile stemness cells, descriptive (cell-level,
        # not a patient-arm group test, so MIN_PATIENTS_PER_ARM does not apply)
        hi_thresh = np.quantile(stemness, 0.9)
        lo_thresh = np.quantile(stemness, 0.1)
        expr_hi = expr[stemness >= hi_thresh]
        expr_lo = expr[stemness <= lo_thresh]
        rs_stat, rs_p = stats.ranksums(expr_hi, expr_lo)

        print(f"    {gene}: Pearson r={pearson_r:.3f} (p={pearson_p:.3g}), "
              f"Spearman rho={spearman_r:.3f} (p={spearman_p:.3g}); "
              f"top-decile-stemness median={np.median(expr_hi):.3f} vs bottom-decile median={np.median(expr_lo):.3f} "
              f"(rank-sum p={rs_p:.3g})")

        rows.append({
            "cohort": name,
            "gene": gene,
            "n_cells": mal.n_obs,
            "n_patients": n_patients,
            "n_stemness_genes_found": len(stem_available),
            "stemness_genes_found": ",".join(stem_available),
            "pearson_r": round(pearson_r, 4),
            "pearson_p": pearson_p,
            "spearman_rho": round(spearman_r, 4),
            "spearman_p": spearman_p,
            "median_top_decile_stemness": round(float(np.median(expr_hi)), 4),
            "median_bottom_decile_stemness": round(float(np.median(expr_lo)), 4),
            "top_vs_bottom_decile_direction": "higher_in_high_stemness" if np.median(expr_hi) > np.median(expr_lo) else "lower_in_high_stemness",
            "top_vs_bottom_decile_p": rs_p,
        })

    del adata, mal
    return pd.DataFrame(rows) if rows else None, plot_data


def make_scatter(name, plot_data):
    """Separate single-gene figures: ACADL -> main Figure 6, CPT1A -> Supplementary
    Figure S2, per the manuscript's decision to report CPT1A only in Supplementary
    Results (pseudoreplication concern on the CPT1A cell-level correlations) rather
    than as a co-equal main-text panel next to ACADL."""
    FIGURE_FILE_STEM = {
        "ACADL": ("Figure_ACADL_stemness_correlation", "ACADL"),
        "CPT1A": ("FigureS2_CPT1A_stemness_correlation", "CPT1A (Supplementary)"),
    }
    for gene, (stem, label) in FIGURE_FILE_STEM.items():
        if gene not in plot_data:
            continue
        fig, ax = plt.subplots(figsize=(5.5, 4.5))
        ax.hexbin(plot_data["stemness_score"], plot_data[gene], gridsize=40, cmap="viridis", mincnt=1)
        ax.set_xlabel("Stemness score (PROM1/CD44/NANOG/KLF4/POU5F1/SOX2)")
        ax.set_ylabel(f"{gene} expression")
        ax.set_title(f"{name}: {label} vs stemness (malignant cells)", fontsize=10)
        plt.tight_layout()
        path = os.path.join(FIGURES_DIR, f"{stem}_{name}.pdf")
        plt.savefig(path, bbox_inches="tight")
        plt.savefig(path.replace(".pdf", ".png"), bbox_inches="tight", dpi=200)
        plt.close()
        print(f"  Scatter saved: {path}")


def main():
    print("=== ACADL/CPT1A vs Stemness Correlation (FAO-stemness literature reconciliation) ===\n")
    cohort_cfg, gene_sets = load_configs()
    stemness_genes = gene_sets["pdac_stemness"]
    print(f"Stemness marker panel (from PMC11351511): {stemness_genes}\n")

    all_results = []
    for cohort in cohort_cfg["cohorts"]:
        print(f"--- {cohort['name']} ---")
        try:
            result_df, plot_data = analyze_cohort(cohort, stemness_genes)
            if result_df is not None:
                all_results.append(result_df)
                make_scatter(cohort["name"], plot_data)
        except Exception as e:
            import traceback
            print(f"  ERROR in {cohort['name']}: {e}")
            traceback.print_exc()
        print()

    if not all_results:
        print("No cohorts produced results. Exiting.")
        return

    combined = pd.concat(all_results, ignore_index=True)
    out_path = os.path.join(TABLES_DIR, "acadl_stemness_correlation.tsv")
    combined.to_csv(out_path, sep="\t", index=False)
    print(f"Results saved: {out_path}\n")

    print("--- Summary ---")
    for gene in FAO_GENES_OF_INTEREST:
        sub = combined[combined["gene"] == gene]
        if sub.empty:
            continue
        n_pos = (sub["pearson_r"] > 0).sum()
        n_sig_pos = ((sub["pearson_r"] > 0) & (sub["pearson_p"] < 0.05)).sum()
        print(f"  {gene}: positive correlation with stemness in {n_pos}/{len(sub)} cohorts "
              f"({n_sig_pos}/{len(sub)} significant, p<0.05)")

    print("\n=== ACADL/CPT1A vs Stemness Correlation Complete ===")


if __name__ == "__main__":
    main()
