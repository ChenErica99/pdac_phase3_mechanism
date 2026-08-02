"""
analyze_acadl_survival.py
ACADL-specific survival analysis (novelty-check recommendation #1, 2026-07-18).

Section 8 of the Phase 3 report tests survival on the composite hypoxia-high/
acinar-low signature, which Section 3 already showed is a poor proxy for a
single coherent cell state (r~0 co-occurrence). Since the manuscript is now
anchored on ACADL specifically (the one protein that independently replicates
FDR<0.05 in both CPTAC centers and survives purity/grade/stage adjustment —
see Section 7), survival should be tested directly against ACADL rather than
against the composite signature it was never meant to summarize.

Cohorts (real data only):
  - TCGA-PAAD  (../Erica (1)/data/processed)                        n=177
  - GSE79668   (../pdac_phase2_validation/data/processed/validation) n=49
  - GSE71729   (../pdac_phase2_validation/data/processed/validation) n=125
  - GSE21501   (data/processed, parsed from SOFT file)               n=102
  - GSE62165 excluded: 0/131 samples have survival_time/event in this
    pipeline (unchanged from Phase 2 — not re-derivable from what was
    downloaded, not a new gap introduced here).

For each cohort:
  1. Continuous Cox (primary): z-scored ACADL as the sole covariate. Reports
     HR per +1 SD of ACADL; HR<1 means higher ACADL is protective, consistent
     with the "ACADL loss -> worse outcome" hypothesis. Avoids the power loss
     of an arbitrary median split.
  2. Median-split Cox + KM (secondary, for visualization and consistency with
     the rest of this pipeline's convention): ACADL-low vs ACADL-high.
Then fixed- and random-effects meta-analysis of the continuous per-SD log(HR)
across the four cohorts.
"""

import os
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

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESEARCH_DIR = os.path.dirname(BASE_DIR)
PHASE2_DIR = os.path.join(RESEARCH_DIR, "pdac_phase2_validation", "data", "processed", "validation")
PHASE1_DIR = os.path.join(RESEARCH_DIR, "pdac_phase1_discovery", "data", "processed")
GSE21501_PARSED = os.path.join(BASE_DIR, "data", "processed", "GSE21501_clinical_expression.tsv")

TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
FIGURES_DIR = os.path.join(BASE_DIR, "results", "figures")
os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def load_gse_cohort(cohort):
    """GSE79668 / GSE71729: raw ACADL expression from the gene-symbol matrix,
    survival from the already-processed sample_scores table (same source
    Phase 2 used for the composite-signature survival test)."""
    scores_path = os.path.join(PHASE2_DIR, f"{cohort}_sample_scores.tsv")
    expr_path = os.path.join(PHASE2_DIR, f"{cohort}_expression_gene_symbol.tsv.gz")
    if not (os.path.exists(scores_path) and os.path.exists(expr_path)):
        print(f"  {cohort}: missing processed files, skipping.")
        return None

    scores = pd.read_csv(scores_path, sep="\t", index_col=0)
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    if "ACADL" not in expr.index:
        print(f"  {cohort}: ACADL not found in expression matrix, skipping.")
        return None

    df = scores[["survival_time_days", "survival_event"]].copy()
    df["ACADL"] = expr.loc["ACADL"].reindex(df.index)
    df = df.dropna(subset=["survival_time_days", "survival_event", "ACADL"])
    df = df[df["survival_time_days"] > 0]
    df.columns = ["T", "E", "ACADL"]
    return df


def load_tcga_paad():
    expr_path = os.path.join(PHASE1_DIR, "tcga_paad_expr_log2_curated.tsv")
    meta_path = os.path.join(PHASE1_DIR, "tcga_paad_sample_metadata_curated.tsv")
    if not (os.path.exists(expr_path) and os.path.exists(meta_path)):
        print("  TCGA-PAAD: missing processed files, skipping.")
        return None

    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)
    if "ACADL" not in expr.index:
        print("  TCGA-PAAD: ACADL not found in expression matrix, skipping.")
        return None

    df = meta[["OS_time_days", "OS_event"]].copy()
    df["ACADL"] = expr.loc["ACADL"].reindex(df.index)
    df = df.dropna(subset=["OS_time_days", "OS_event", "ACADL"])
    df = df[df["OS_time_days"] > 0]
    df.columns = ["T", "E", "ACADL"]
    return df


def load_gse21501():
    if not os.path.exists(GSE21501_PARSED):
        print("  GSE21501: parsed table not found — run parse_gse21501_survival.py first.")
        return None
    df = pd.read_csv(GSE21501_PARSED, sep="\t")
    df["os_time"] = pd.to_numeric(df["os_time"], errors="coerce")
    df["os_event"] = pd.to_numeric(df["os_event"], errors="coerce")
    df = df.dropna(subset=["os_time", "os_event", "ACADL"])
    df = df[df["os_time"] > 0]
    df = df[["os_time", "os_event", "ACADL"]]
    df.columns = ["T", "E", "ACADL"]
    return df


def check_ph_assumption(cph, fit_df, covariate):
    """Schoenfeld-residual proportional-hazards check (Limitations item: not
    yet performed for any Cox model, flagged 2026-07-25/26). Returns the
    rank-transform test p-value for the covariate (lifelines default
    transform) plus a global (any-transform) verdict; a low p-value (<0.05)
    indicates a covariate effect that changes over time, i.e. a PH violation."""
    from lifelines.statistics import proportional_hazard_test
    result = proportional_hazard_test(cph, fit_df, time_transform="rank")
    row = result.summary.loc[covariate]
    return {
        "ph_test_p": round(float(row["p"]), 5),
        "ph_test_stat": round(float(row["test_statistic"]), 4),
        "ph_violation": bool(row["p"] < 0.05),
    }


def analyze_cohort(name, df):
    """Fit continuous and median-split Cox models; return a result dict plus
    the annotated dataframe (for the KM plot)."""
    from lifelines import CoxPHFitter
    from lifelines.statistics import logrank_test

    n_total = len(df)
    n_events = int(df["E"].sum())
    df = df.copy()
    df["ACADL_z"] = (df["ACADL"] - df["ACADL"].mean()) / df["ACADL"].std()
    med = df["ACADL"].median()
    df["ACADL_low"] = (df["ACADL"] < med).astype(int)
    n_low = int(df["ACADL_low"].sum())
    print(f"  {name}: n={n_total}, events={n_events}, ACADL-low={n_low}, ACADL-high={n_total - n_low}")

    if n_total < 10 or n_events < 5:
        print(f"  {name}: insufficient events for a Cox fit, skipping.")
        return None, df

    # 1. Continuous Cox (primary)
    cph_c = CoxPHFitter()
    cph_c.fit(df[["T", "E", "ACADL_z"]], duration_col="T", event_col="E")
    s_c = cph_c.summary
    hr_c = float(s_c.loc["ACADL_z", "exp(coef)"])
    hr_c_lo = float(s_c.loc["ACADL_z", "exp(coef) lower 95%"])
    hr_c_hi = float(s_c.loc["ACADL_z", "exp(coef) upper 95%"])
    p_c = float(s_c.loc["ACADL_z", "p"])
    log_hr_c = float(s_c.loc["ACADL_z", "coef"])
    se_c = float(s_c.loc["ACADL_z", "se(coef)"])
    ph_c = check_ph_assumption(cph_c, df[["T", "E", "ACADL_z"]], "ACADL_z")

    # 2. Median-split Cox + log-rank (secondary)
    cph_m = CoxPHFitter()
    cph_m.fit(df[["T", "E", "ACADL_low"]], duration_col="T", event_col="E")
    s_m = cph_m.summary
    hr_m = float(s_m.loc["ACADL_low", "exp(coef)"])
    hr_m_lo = float(s_m.loc["ACADL_low", "exp(coef) lower 95%"])
    hr_m_hi = float(s_m.loc["ACADL_low", "exp(coef) upper 95%"])
    p_m = float(s_m.loc["ACADL_low", "p"])
    ph_m = check_ph_assumption(cph_m, df[["T", "E", "ACADL_low"]], "ACADL_low")

    lr = logrank_test(df.loc[df["ACADL_low"] == 1, "T"], df.loc[df["ACADL_low"] == 0, "T"],
                       event_observed_A=df.loc[df["ACADL_low"] == 1, "E"],
                       event_observed_B=df.loc[df["ACADL_low"] == 0, "E"])

    print(f"    Continuous:   HR/SD={hr_c:.3f} [{hr_c_lo:.3f}-{hr_c_hi:.3f}] p={p_c:.4f}  "
          f"PH test p={ph_c['ph_test_p']:.4f}{' (VIOLATION)' if ph_c['ph_violation'] else ''}")
    print(f"    Median-split: HR(low vs high)={hr_m:.3f} [{hr_m_lo:.3f}-{hr_m_hi:.3f}] p={p_m:.4f} (log-rank p={lr.p_value:.4f})  "
          f"PH test p={ph_m['ph_test_p']:.4f}{' (VIOLATION)' if ph_m['ph_violation'] else ''}")

    return {
        "cohort": name,
        "n_total": n_total,
        "n_events": n_events,
        "n_ACADL_low": n_low,
        "hr_per_sd": round(hr_c, 4),
        "hr_per_sd_lower95": round(hr_c_lo, 4),
        "hr_per_sd_upper95": round(hr_c_hi, 4),
        "p_continuous": round(p_c, 5),
        "log_hr_per_sd": log_hr_c,
        "se_log_hr_per_sd": se_c,
        "hr_low_vs_high": round(hr_m, 4),
        "hr_low_vs_high_lower95": round(hr_m_lo, 4),
        "hr_low_vs_high_upper95": round(hr_m_hi, 4),
        "p_median_split": round(p_m, 5),
        "logrank_p": round(lr.p_value, 5),
        "ph_test_p_continuous": ph_c["ph_test_p"],
        "ph_violation_continuous": ph_c["ph_violation"],
        "ph_test_p_median_split": ph_m["ph_test_p"],
        "ph_violation_median_split": ph_m["ph_violation"],
        "is_simulated": False,
    }, df


def make_km_plot(name, df):
    from lifelines import KaplanMeierFitter
    fig, ax = plt.subplots(figsize=(6, 4.5))
    kmf_lo = KaplanMeierFitter()
    kmf_hi = KaplanMeierFitter()
    low = df[df["ACADL_low"] == 1]
    high = df[df["ACADL_low"] == 0]
    kmf_lo.fit(low["T"], event_observed=low["E"], label=f"ACADL-low (n={len(low)})")
    kmf_hi.fit(high["T"], event_observed=high["E"], label=f"ACADL-high (n={len(high)})")
    kmf_lo.plot_survival_function(ax=ax, color="#c0392b", ci_show=True)
    kmf_hi.plot_survival_function(ax=ax, color="#2980b9", ci_show=True)
    ax.set_title(f"{name}: overall survival by ACADL level", fontsize=10)
    ax.set_xlabel("Time")
    ax.set_ylabel("Survival probability")
    ax.legend(fontsize=8)
    plt.tight_layout()
    path = os.path.join(FIGURES_DIR, f"KM_{name}_ACADL_low_vs_high.pdf")
    plt.savefig(path, bbox_inches="tight")
    plt.close()
    print(f"  KM plot saved: {path}")


def fixed_effects_meta(log_hrs, ses):
    weights = 1.0 / (np.array(ses) ** 2)
    pooled_log_hr = np.sum(weights * np.array(log_hrs)) / np.sum(weights)
    pooled_se = np.sqrt(1.0 / np.sum(weights))
    z = pooled_log_hr / pooled_se
    p_val = 2 * (1 - stats.norm.cdf(abs(z)))
    hr = np.exp(pooled_log_hr)
    return {
        "model": "fixed_effects", "HR_pooled": round(hr, 4),
        "HR_lower_95CI": round(np.exp(pooled_log_hr - 1.96 * pooled_se), 4),
        "HR_upper_95CI": round(np.exp(pooled_log_hr + 1.96 * pooled_se), 4),
        "p_value": round(p_val, 5),
    }


def random_effects_meta(log_hrs, ses):
    log_hrs = np.array(log_hrs)
    ses = np.array(ses)
    k = len(log_hrs)
    weights_fe = 1.0 / ses ** 2
    mu_fe = np.sum(weights_fe * log_hrs) / np.sum(weights_fe)
    Q = np.sum(weights_fe * (log_hrs - mu_fe) ** 2)
    df_ = k - 1
    C = np.sum(weights_fe) - np.sum(weights_fe ** 2) / np.sum(weights_fe)
    tau2 = max(0.0, (Q - df_) / C) if C > 0 else 0.0
    weights_re = 1.0 / (ses ** 2 + tau2)
    mu_re = np.sum(weights_re * log_hrs) / np.sum(weights_re)
    se_re = np.sqrt(1.0 / np.sum(weights_re))
    z = mu_re / se_re
    p_val = 2 * (1 - stats.norm.cdf(abs(z)))
    I2 = max(0.0, (Q - df_) / Q * 100) if Q > 0 and df_ > 0 else 0.0
    p_het = 1 - stats.chi2.cdf(Q, df_) if df_ > 0 else np.nan
    return {
        "model": "random_effects_DL", "HR_pooled": round(np.exp(mu_re), 4),
        "HR_lower_95CI": round(np.exp(mu_re - 1.96 * se_re), 4),
        "HR_upper_95CI": round(np.exp(mu_re + 1.96 * se_re), 4),
        "p_value": round(p_val, 5), "tau2": round(tau2, 4),
        "I2_pct": round(I2, 1), "p_heterogeneity": round(p_het, 4) if not np.isnan(p_het) else np.nan,
    }


def make_forest_plot(results, meta_fe, meta_re):
    n = len(results)
    fig, ax = plt.subplots(figsize=(9, max(4, n * 0.9 + 3)))
    y_positions = list(range(n, 0, -1))
    for r, y in zip(results, y_positions):
        ax.plot([r["hr_per_sd_lower95"], r["hr_per_sd_upper95"]], [y, y], color="#2c3e50", linewidth=1.5)
        ax.plot(r["hr_per_sd"], y, "s", color="#2c3e50",
                 markersize=8 + 4 * np.sqrt(r["n_events"] / 80))
        ax.text(-0.1, y, r["cohort"], ha="right", va="center", fontsize=9, transform=ax.get_yaxis_transform())
        ax.text(1.02, y, f'HR/SD={r["hr_per_sd"]:.2f} [{r["hr_per_sd_lower95"]:.2f}-{r["hr_per_sd_upper95"]:.2f}] p={r["p_continuous"]:.3f}',
                ha="left", va="center", fontsize=8, transform=ax.get_yaxis_transform())

    for y_m, meta, label, color in [(0.4, meta_fe, "Pooled (FE)", "#e74c3c"), (-0.1, meta_re, "Pooled (RE)", "#8e44ad")]:
        ax.plot([meta["HR_lower_95CI"], meta["HR_upper_95CI"]], [y_m, y_m], color=color, linewidth=2)
        ax.plot(meta["HR_pooled"], y_m, "D", color=color, markersize=10)
        ax.text(-0.1, y_m, label, ha="right", va="center", fontsize=9, color=color, transform=ax.get_yaxis_transform())
        ax.text(1.02, y_m, f'HR/SD={meta["HR_pooled"]:.2f} [{meta["HR_lower_95CI"]:.2f}-{meta["HR_upper_95CI"]:.2f}] p={meta["p_value"]:.3f}',
                ha="left", va="center", fontsize=8, color=color, transform=ax.get_yaxis_transform())

    ax.axvline(1.0, color="black", linestyle="--", linewidth=1, alpha=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("Hazard Ratio per +1 SD of ACADL (log scale) — HR<1 = higher ACADL protective", fontsize=9)
    ax.set_xlim(0.3, 3.0)
    ax.set_ylim(-0.5, n + 1)
    ax.set_yticks([])
    ax.set_title("ACADL-Specific Survival: Continuous Cox HR per Cohort\n"
                 f"I²={meta_re.get('I2_pct', 0):.0f}%, p_het={meta_re.get('p_heterogeneity', 1):.3f}", fontsize=11)
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    path = os.path.join(FIGURES_DIR, "Figure_ACADL_survival_forest.pdf")
    plt.savefig(path, bbox_inches="tight")
    plt.close()
    print(f"\nForest plot saved: {path}")


def main():
    print("=== ACADL-Specific Survival Analysis ===\n")

    loaders = {
        "TCGA-PAAD": load_tcga_paad,
        "GSE79668": lambda: load_gse_cohort("GSE79668"),
        "GSE71729": lambda: load_gse_cohort("GSE71729"),
        "GSE21501": load_gse21501,
    }

    results = []
    for name, loader in loaders.items():
        print(f"\n--- {name} ---")
        df = loader()
        if df is None or len(df) == 0:
            print(f"  {name}: no usable data.")
            continue
        result, annotated_df = analyze_cohort(name, df)
        if result is not None:
            results.append(result)
            make_km_plot(name, annotated_df)

    print(f"\n--- GSE62165 ---\n  Excluded: 0/131 samples have survival_time/event in this pipeline "
          f"(unchanged since Phase 2 — the array series has no survival annotation available).")

    if results:
        print("\n--- Proportional-hazards assumption (Schoenfeld residuals, rank transform) ---")
        n_viol_c = sum(r["ph_violation_continuous"] for r in results)
        n_viol_m = sum(r["ph_violation_median_split"] for r in results)
        for r in results:
            flag_c = " <-- VIOLATION" if r["ph_violation_continuous"] else ""
            flag_m = " <-- VIOLATION" if r["ph_violation_median_split"] else ""
            print(f"  {r['cohort']}: continuous p={r['ph_test_p_continuous']:.4f}{flag_c}; "
                  f"median-split p={r['ph_test_p_median_split']:.4f}{flag_m}")
        print(f"  Summary: {n_viol_c}/{len(results)} continuous models and {n_viol_m}/{len(results)} "
              f"median-split models show a PH violation at p<0.05.")

    if len(results) == 0:
        print("\nNo cohorts had sufficient data for a Cox fit. Exiting.")
        return

    results_df = pd.DataFrame(results)
    out_path = os.path.join(TABLES_DIR, "acadl_survival_by_cohort.tsv")
    results_df.to_csv(out_path, sep="\t", index=False)
    print(f"\nPer-cohort results saved: {out_path}")

    if len(results) >= 2:
        log_hrs = [r["log_hr_per_sd"] for r in results]
        ses = [r["se_log_hr_per_sd"] for r in results]
        meta_fe = fixed_effects_meta(log_hrs, ses)
        meta_re = random_effects_meta(log_hrs, ses)

        print("\nFixed-effects meta (continuous, per +1 SD ACADL):")
        print(f"  HR={meta_fe['HR_pooled']:.3f} [{meta_fe['HR_lower_95CI']:.3f}-{meta_fe['HR_upper_95CI']:.3f}] p={meta_fe['p_value']:.4f}")
        print("Random-effects meta:")
        print(f"  HR={meta_re['HR_pooled']:.3f} [{meta_re['HR_lower_95CI']:.3f}-{meta_re['HR_upper_95CI']:.3f}] "
              f"p={meta_re['p_value']:.4f}, I²={meta_re.get('I2_pct', 0):.0f}%")

        meta_rows = []
        for meta, label in [(meta_fe, "pooled_fixed_effects"), (meta_re, "pooled_random_effects")]:
            row = {"model": label, "n_cohorts": len(results), **meta}
            meta_rows.append(row)
        meta_df = pd.DataFrame(meta_rows)
        meta_out = os.path.join(TABLES_DIR, "acadl_survival_meta_analysis.tsv")
        meta_df.to_csv(meta_out, sep="\t", index=False)
        print(f"Meta-analysis saved: {meta_out}")

        make_forest_plot(results, meta_fe, meta_re)
    else:
        print("\nFewer than 2 cohorts with valid fits — skipping meta-analysis.")

    print("\n=== ACADL Survival Analysis Complete ===")


if __name__ == "__main__":
    main()
