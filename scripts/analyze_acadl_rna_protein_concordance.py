"""
analyze_acadl_rna_protein_concordance.py
RNA-protein concordance for ACADL specifically (novelty-check recommendation
#3, 2026-07-18).

Section 7 established that ACADL protein is reduced in aggressive tumors,
replicated FDR<0.05 in both CPTAC centers, and survives purity/grade/stage
adjustment. It never checked whether that reduction is also visible at the
mRNA level in CPTAC's own matched washu transcriptomics (already loaded for
group assignment, so this is essentially free). Two questions:

  1. Do ACADL protein and ACADL mRNA correlate across tumor samples at all
     (person-level concordance, not just group-level)? A weak/absent
     correlation would point to post-transcriptional regulation of ACADL
     abundance, distinguishing it mechanistically from the SREBP lipogenic
     genes (Section 7 already noted SREBF1 itself isn't elevated at protein
     level, suggesting post-translational SREBP1 activity regulation there).
  2. Is ACADL mRNA itself reduced in aggressive vs reference tumors (same
     group definition Section 7 uses, which is built from hypoxia/acinar
     genes, not ACADL, so this test is not circular)?

Reuses the exact CPTAC loading/group-assignment code from
analyze_cptac_protein.py rather than re-implementing it, per this pipeline's
"run the real production function, don't reimplement" lesson (see Phase 3
memory on the CAF subtype groupby bug).
"""

import os
import sys
import numpy as np
import pandas as pd
import random
import warnings
warnings.filterwarnings("ignore")

np.random.seed(1234)
random.seed(1234)

from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_cptac_protein import load_cptac_dataset, load_proteomics, assign_groups

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
FIGURES_DIR = os.path.join(BASE_DIR, "results", "figures")
os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

GENE = "ACADL"


def get_gene_series(df, gene):
    """CPTAC transcriptomics has MultiIndex columns (Name, Database_ID);
    df[gene] returns a 1-column DataFrame, not a Series, in that case."""
    col = df[gene]
    if isinstance(col, pd.DataFrame):
        col = col.mean(axis=1)
    return col


def concordance_for_source(source, prot, trans_acadl, group_df):
    common = prot.index.intersection(trans_acadl.index).intersection(group_df.index)
    print(f"  {source}: {len(common)} samples with protein + RNA + group assignment")

    prot_g = prot.loc[common, GENE]
    trans_g = trans_acadl.loc[common]
    both = pd.DataFrame({"protein": prot_g, "rna": trans_g}).dropna()

    rho, p_corr = stats.spearmanr(both["protein"], both["rna"])
    print(f"    RNA-protein Spearman rho={rho:.3f}, p={p_corr:.4g}, n={len(both)}")

    grp = group_df.loc[common, "group"]
    agg_idx = common[grp == "aggressive"]
    ref_idx = common[grp == "reference"]

    prot_agg = prot.loc[agg_idx, GENE].dropna()
    prot_ref = prot.loc[ref_idx, GENE].dropna()
    stat_p, pval_p = stats.ranksums(prot_agg, prot_ref)
    dir_p = "down" if prot_agg.median() < prot_ref.median() else "up"

    rna_agg = trans_acadl.loc[agg_idx].dropna()
    rna_ref = trans_acadl.loc[ref_idx].dropna()
    stat_r, pval_r = stats.ranksums(rna_agg, rna_ref)
    dir_r = "down" if rna_agg.median() < rna_ref.median() else "up"

    print(f"    Protein aggressive vs reference: {dir_p}, p={pval_p:.4g} (n={len(prot_agg)} vs {len(prot_ref)})")
    print(f"    RNA aggressive vs reference:     {dir_r}, p={pval_r:.4g} (n={len(rna_agg)} vs {len(rna_ref)})")

    # Scatter: protein vs RNA colored by group
    fig, ax = plt.subplots(figsize=(5.5, 5))
    plot_df = both.join(grp.rename("group"))
    palette = {"aggressive": "#e74c3c", "reference": "#3498db", "other": "#bdc3c7"}
    for g, color in palette.items():
        sub = plot_df[plot_df["group"] == g]
        ax.scatter(sub["rna"], sub["protein"], color=color, alpha=0.6, s=18, label=f"{g} (n={len(sub)})")
    ax.set_xlabel(f"{GENE} mRNA (washu transcriptomics)")
    ax.set_ylabel(f"{GENE} protein ({source})")
    ax.set_title(f"{GENE} RNA-protein concordance — {source}\nSpearman rho={rho:.3f}, p={p_corr:.3g}", fontsize=10)
    ax.legend(fontsize=7)
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, f"Figure_ACADL_RNA_protein_scatter_{source}.pdf")
    plt.savefig(fig_path, bbox_inches="tight")
    plt.savefig(fig_path.replace(".pdf", ".png"), bbox_inches="tight", dpi=200)
    plt.close()
    print(f"    Scatter saved: {fig_path}")

    return {
        "source": source,
        "gene": GENE,
        "n_paired_rna_protein": len(both),
        "rna_protein_spearman_rho": round(rho, 4),
        "rna_protein_spearman_p": pval_from_round(p_corr),
        "protein_n_aggressive": len(prot_agg),
        "protein_n_reference": len(prot_ref),
        "protein_median_aggressive": round(float(prot_agg.median()), 4),
        "protein_median_reference": round(float(prot_ref.median()), 4),
        "protein_direction": dir_p,
        "protein_p_value": pval_from_round(pval_p),
        "rna_n_aggressive": len(rna_agg),
        "rna_n_reference": len(rna_ref),
        "rna_median_aggressive": round(float(rna_agg.median()), 4),
        "rna_median_reference": round(float(rna_ref.median()), 4),
        "rna_direction": dir_r,
        "rna_p_value": pval_from_round(pval_r),
        "rna_and_protein_concordant_direction": dir_p == dir_r,
    }


def pval_from_round(p):
    return round(float(p), 6)


def main():
    print("=== ACADL RNA-Protein Concordance (CPTAC-PDA) ===\n")
    pdac, trans, purity_df, clin = load_cptac_dataset()

    if GENE not in trans.columns:
        print(f"{GENE} not found in washu transcriptomics — cannot proceed.")
        return

    print("\nAssigning groups from transcriptomics (same definition as Section 7)...")
    group_df = assign_groups(trans)
    trans_acadl = get_gene_series(trans, GENE)

    results = []
    for source in ["umich", "bcm"]:
        print(f"\n--- {source} ---")
        prot = load_proteomics(pdac, source)
        if GENE not in prot.columns:
            print(f"  {GENE} not found in {source} proteomics, skipping.")
            continue
        results.append(concordance_for_source(source, prot, trans_acadl, group_df))

    if not results:
        print("\nNo sources had ACADL protein data. Exiting.")
        return

    results_df = pd.DataFrame(results)
    out_path = os.path.join(TABLES_DIR, "acadl_rna_protein_concordance.tsv")
    results_df.to_csv(out_path, sep="\t", index=False)
    print(f"\nResults saved: {out_path}")

    print("\n--- Summary ---")
    for _, row in results_df.iterrows():
        print(f"  {row['source']}: RNA-protein rho={row['rna_protein_spearman_rho']:.3f} "
              f"(p={row['rna_protein_spearman_p']:.3g}); "
              f"RNA {row['rna_direction']} (p={row['rna_p_value']:.3g}) vs "
              f"protein {row['protein_direction']} (p={row['protein_p_value']:.3g}) "
              f"aggressive-vs-reference — concordant: {row['rna_and_protein_concordant_direction']}")

    print("\n=== ACADL RNA-Protein Concordance Complete ===")


if __name__ == "__main__":
    main()
