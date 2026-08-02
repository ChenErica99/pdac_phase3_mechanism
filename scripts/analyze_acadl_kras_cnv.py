"""
analyze_acadl_kras_cnv.py
Reviewer-anticipated analyses (2026-08-02): the manuscript's Discussion cites
Nwosu et al. 2025 (ACADL downregulation reported in a KRAS-mutant PDAC
context) and lists ACADL-locus copy-number characterization as a "future
direction." Both are directly testable in data already loaded elsewhere in
this pipeline (cptac washu somatic_mutation and washu CNV), so this script
runs them now rather than deferring them:

  1. KRAS-mutation stratification: does ACADL suppression in the aggressive
     group hold within the KRAS-mutant subset specifically? (Nwosu et al.'s
     finding was in a KRAS-mutant context; almost all CPTAC-PDA patients are
     KRAS-mutant, so this mainly checks that the effect isn't somehow
     confined to the small KRAS-WT minority.)
  2. ACADL copy number: is ACADL protein/RNA suppression accompanied by
     copy-number loss at the ACADL locus, or is it purely
     transcriptional/epigenetic? Tests CNV vs protein/RNA correlation, CNV
     by aggressive/reference group, and CNV as an added covariate in the
     continuous confound model (does it explain away the acinar effect?).

Reuses load_cptac_dataset/load_proteomics/assign_groups/zscore_genes from
analyze_cptac_protein.py rather than reimplementing CPTAC loading.
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
import statsmodels.api as sm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_cptac_protein import (
    load_cptac_dataset, load_proteomics, assign_groups, zscore_genes,
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


def load_kras_status(pdac):
    """Binary KRAS-mutant status per patient from the washu somatic MAF."""
    print("  Fetching washu somatic mutation data...")
    mut = pdac.get_somatic_mutation(source="washu")
    kras_patients = set(mut.loc[mut["Gene"] == "KRAS"].index.unique())
    all_patients = mut.index.unique()
    status = pd.Series(
        [p in kras_patients for p in all_patients], index=all_patients, name="kras_mutant"
    )
    print(f"  KRAS-mutant: {status.sum()}/{len(status)} patients "
          f"({100 * status.mean():.1f}%)")
    return status


def load_cnv(pdac, source, gene):
    """Gene-level log2 copy-number ratio for one gene, one source."""
    print(f"  Fetching {source} CNV...")
    cnv_raw = pdac.get_CNV(source=source)
    cnv = cnv_raw[~cnv_raw.index.str.endswith(".N")].copy()
    if gene not in cnv.columns:
        print(f"  {gene} not found in {source} CNV.")
        return None
    return get_gene_series(cnv, gene)


def kras_stratified_test(group_df, prot, rna, source, kras_status):
    """Aggressive-vs-reference test for ACADL, restricted to KRAS-mutant patients."""
    common = group_df.index.intersection(kras_status.index)
    mutant_idx = common[kras_status.loc[common]]
    wt_idx = common[~kras_status.loc[common]]

    rows = []
    for label, values, val_name in [("protein", prot, source), ("rna", rna, "washu_rna")]:
        if values is None:
            continue
        agg_all = group_df.index[group_df["group"] == "aggressive"]
        ref_all = group_df.index[group_df["group"] == "reference"]

        for subset_name, subset_idx in [("kras_mutant", mutant_idx), ("kras_wt", wt_idx)]:
            agg = values.reindex(agg_all).loc[values.reindex(agg_all).index.intersection(subset_idx)].dropna()
            ref = values.reindex(ref_all).loc[values.reindex(ref_all).index.intersection(subset_idx)].dropna()
            if len(agg) < 3 or len(ref) < 3:
                rows.append({
                    "source": val_name, "subset": subset_name, "n_aggressive": len(agg),
                    "n_reference": len(ref), "median_aggressive": np.nan,
                    "median_reference": np.nan, "direction": None, "p_value": None,
                    "note": "too few samples to test",
                })
                continue
            stat, p = stats.ranksums(agg, ref)
            direction = "down" if agg.median() < ref.median() else "up"
            rows.append({
                "source": val_name, "subset": subset_name, "n_aggressive": len(agg),
                "n_reference": len(ref), "median_aggressive": round(float(agg.median()), 4),
                "median_reference": round(float(ref.median()), 4), "direction": direction,
                "p_value": round(float(p), 6), "note": "",
            })
    return rows


def cnv_confound_model(y, cnv, scores, purity_df, clin, label):
    common = y.index.intersection(cnv.index).intersection(scores.index)
    if purity_df is not None:
        common = common.intersection(purity_df.index)
    if clin is not None:
        common = common.intersection(clin.index)

    df = pd.DataFrame(index=common)
    df["y"] = y.loc[common].astype(float)
    df["cnv"] = cnv.loc[common].astype(float)
    df["hypoxia_z"] = scores.loc[common, "hypoxia_z"]
    df["acinar_z"] = scores.loc[common, "acinar_z"]
    df["purity"] = purity_df.loc[common, "TumorPurity"] if purity_df is not None else np.nan
    df["grade"] = clin.loc[common, "grade_numeric"] if clin is not None else np.nan
    df["stage"] = clin.loc[common, "stage_numeric"] if clin is not None else np.nan

    df = df.dropna(subset=["y", "cnv"])
    for c in ["purity", "grade", "stage"]:
        if df[c].notna().any():
            df[c] = df[c].fillna(df[c].median())
    df = df.dropna()

    n = len(df)
    print(f"  [{label}] n={n} samples with complete data (incl. CNV)")
    if n < 15:
        print(f"  [{label}] too few samples, skipping.")
        return None

    # CNV-only correlation
    rho, p_rho = stats.spearmanr(df["y"], df["cnv"])

    X = sm.add_constant(df[["hypoxia_z", "acinar_z", "purity", "grade", "stage", "cnv"]])
    model = sm.OLS(df["y"], X).fit()

    row = {"label": label, "n": n, "spearman_rho_y_vs_cnv": round(rho, 4),
           "spearman_p_y_vs_cnv": round(float(p_rho), 6)}
    for term in ["hypoxia_z", "acinar_z", "purity", "grade", "stage", "cnv"]:
        row[f"coef_{term}"] = round(float(model.params[term]), 5)
        row[f"pval_{term}"] = round(float(model.pvalues[term]), 6)
    row["r_squared"] = round(float(model.rsquared), 4)
    return row


def main():
    print("=== ACADL: KRAS-Mutation Stratification and Copy-Number Analysis ===\n")
    pdac, trans, purity_df, clin = load_cptac_dataset()

    if GENE not in trans.columns:
        print(f"{GENE} not in washu transcriptomics, aborting.")
        return

    print("\nAssigning groups from transcriptomics (same definition as Section 7)...")
    group_df = assign_groups(trans)
    rna = get_gene_series(trans, GENE)

    scores = pd.DataFrame({
        "hypoxia_z": zscore_genes(trans[[g for g in HYPOXIA_GENES if g in trans.columns]]).mean(axis=1),
        "acinar_z": zscore_genes(trans[[g for g in ACINAR_GENES if g in trans.columns]]).mean(axis=1),
    })

    # --- 1. KRAS mutation status ---
    print("\n--- KRAS mutation status ---")
    kras_status = load_kras_status(pdac)

    # Sanity check: does KRAS-mutant rate differ between aggressive/reference groups?
    common = group_df.index.intersection(kras_status.index)
    agg_idx = common[group_df.loc[common, "group"] == "aggressive"]
    ref_idx = common[group_df.loc[common, "group"] == "reference"]
    agg_kras_rate = kras_status.loc[agg_idx].mean()
    ref_kras_rate = kras_status.loc[ref_idx].mean()
    print(f"  KRAS-mutant rate: aggressive={agg_kras_rate:.3f} (n={len(agg_idx)}), "
          f"reference={ref_kras_rate:.3f} (n={len(ref_idx)})")
    table_2x2 = [
        [kras_status.loc[agg_idx].sum(), (~kras_status.loc[agg_idx]).sum()],
        [kras_status.loc[ref_idx].sum(), (~kras_status.loc[ref_idx]).sum()],
    ]
    odds_ratio, fisher_p = stats.fisher_exact(table_2x2)
    print(f"  Fisher's exact (aggressive vs reference KRAS-mutant rate): OR={odds_ratio:.3f}, p={fisher_p:.4g}")

    kras_rate_row = pd.DataFrame([{
        "n_total": len(common), "n_kras_mutant": int(kras_status.loc[common].sum()),
        "pct_kras_mutant": round(100 * kras_status.loc[common].mean(), 1),
        "agg_kras_rate": round(float(agg_kras_rate), 4), "ref_kras_rate": round(float(ref_kras_rate), 4),
        "fisher_or": round(float(odds_ratio), 4), "fisher_p": round(float(fisher_p), 6),
    }])
    kras_rate_row.to_csv(os.path.join(TABLES_DIR, "kras_mutation_rate_by_group.tsv"), sep="\t", index=False)

    print("\n--- ACADL protein/RNA, aggressive vs reference, stratified by KRAS status ---")
    all_rows = []
    for source in ["umich", "bcm"]:
        prot = load_proteomics(pdac, source)
        prot_g = get_gene_series(prot, GENE) if GENE in prot.columns else None
        rows = kras_stratified_test(group_df, prot_g, None, source, kras_status)
        all_rows.extend(rows)
    rna_rows = kras_stratified_test(group_df, None, rna, "washu_rna", kras_status)
    all_rows.extend(rna_rows)
    kras_strat_df = pd.DataFrame(all_rows)
    kras_strat_path = os.path.join(TABLES_DIR, "acadl_kras_stratified.tsv")
    kras_strat_df.to_csv(kras_strat_path, sep="\t", index=False)
    print(kras_strat_df.to_string(index=False))
    print(f"  Saved: {kras_strat_path}")

    # --- 2. ACADL copy number ---
    print("\n--- ACADL copy number (washu, bcm) ---")
    cnv_rows_group = []
    cnv_rows_model = []
    for cnv_source in ["washu", "bcm"]:
        try:
            cnv = load_cnv(pdac, cnv_source, GENE)
        except Exception as e:
            print(f"  {cnv_source} CNV load failed: {e}")
            continue
        if cnv is None:
            continue

        # CNV by aggressive/reference group
        common_g = group_df.index.intersection(cnv.index)
        agg_cnv = cnv.loc[common_g[group_df.loc[common_g, "group"] == "aggressive"]].dropna()
        ref_cnv = cnv.loc[common_g[group_df.loc[common_g, "group"] == "reference"]].dropna()
        stat, p = stats.ranksums(agg_cnv, ref_cnv)
        direction = "lower_(loss)" if agg_cnv.median() < ref_cnv.median() else "higher"
        print(f"  [{cnv_source}] ACADL CNV aggressive (n={len(agg_cnv)}, median={agg_cnv.median():.4f}) "
              f"vs reference (n={len(ref_cnv)}, median={ref_cnv.median():.4f}): {direction}, p={p:.4g}")
        cnv_rows_group.append({
            "cnv_source": cnv_source, "n_aggressive": len(agg_cnv), "n_reference": len(ref_cnv),
            "median_cnv_aggressive": round(float(agg_cnv.median()), 5),
            "median_cnv_reference": round(float(ref_cnv.median()), 5),
            "direction_in_aggressive": direction, "p_value": round(float(p), 6),
        })

        # CNV vs protein/RNA + confound model (washu CNV matched to washu RNA and both proteomics;
        # bcm CNV matched to bcm protein only)
        if cnv_source == "washu":
            targets = [("washu_rna", rna)]
            for prot_source in ["umich", "bcm"]:
                prot = load_proteomics(pdac, prot_source)
                if GENE in prot.columns:
                    targets.append((f"{prot_source}_protein", get_gene_series(prot, GENE)))
        else:
            prot = load_proteomics(pdac, cnv_source)
            targets = [(f"{cnv_source}_protein", get_gene_series(prot, GENE))] if GENE in prot.columns else []

        for label, y in targets:
            print(f"\n  Confound model with CNV covariate: {label} ~ hypoxia_z + acinar_z + purity + grade + stage + cnv({cnv_source})")
            row = cnv_confound_model(y, cnv, scores, purity_df, clin, f"{label}_vs_{cnv_source}_cnv")
            if row is not None:
                print(f"    n={row['n']}  CNV-y Spearman rho={row['spearman_rho_y_vs_cnv']:.3f} (p={row['spearman_p_y_vs_cnv']:.4g})  "
                      f"coef_cnv={row['coef_cnv']:.4f} (p={row['pval_cnv']:.4g})  "
                      f"coef_acinar_z={row['coef_acinar_z']:.4f} (p={row['pval_acinar_z']:.4g})")
                cnv_rows_model.append(row)

    pd.DataFrame(cnv_rows_group).to_csv(os.path.join(TABLES_DIR, "acadl_cnv_by_group.tsv"), sep="\t", index=False)
    pd.DataFrame(cnv_rows_model).to_csv(os.path.join(TABLES_DIR, "acadl_cnv_confound_model.tsv"), sep="\t", index=False)

    print("\n=== Complete ===")
    print("Outputs: results/tables/kras_mutation_rate_by_group.tsv, acadl_kras_stratified.tsv, "
          "acadl_cnv_by_group.tsv, acadl_cnv_confound_model.tsv")


if __name__ == "__main__":
    main()
