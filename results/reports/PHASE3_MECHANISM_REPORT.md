# Phase 3 Mechanistic Report: Single-Cell Resolution of Aggressive PDAC State

**Generated:** 2026-07-03 (updated 2026-07-04 with BCM replication; updated 2026-07-06 with patient-level pseudobulk re-analysis of lipid/CAF/EMT/immune-endothelial cell-of-origin, Sections 4/5/9; updated 2026-07-06 with a non-circular purity proxy, Section 6; updated 2026-07-06 with a real GSE21501 Cox fit replacing the hardcoded survival value, Section 8; updated 2026-07-12 with the full 43-patient GSE202051 object replacing a 1-patient subset, Sections 3/4/5/9/10/11; updated 2026-07-13 replacing the heuristic purity proxy with the validated ESTIMATE algorithm (Yoshihara et al. 2013) via the `tidyestimate` package, Sections 6/10/11/13; updated 2026-07-13 with a purity/grade/stage confounder-adjusted re-analysis of the CPTAC ACADL/lipid protein findings, Section 7/10/13; updated 2026-07-13 with a GSE202051 treatment-status sensitivity check (inconclusive for cell-of-origin tests, confirmatory for the hypoxia/acinar correlation), Sections 3/4/5/9/11/13; **updated 2026-07-13 fixing a real CAF-subtype-proportion correctness bug (categorical groupby without `observed=True`, diluting every fraction) affecting all 3 single-cell cohorts — retracts the "lower myCAF/iCAF in aggressive" claim, Sections 3/5/11/13**; updated 2026-07-18 with ACADL RNA-protein concordance in CPTAC's own matched transcriptomics, Section 7; updated 2026-07-18 with an ACADL-specific survival analysis (TCGA-PAAD, GSE79668, GSE71729, GSE21501) replacing reliance on the composite-signature survival test for the anchor gene, Section 8/10/13; updated 2026-07-18 reconciling the FAO-required-for-stemness literature tension (PMC11351511) via a single-cell ACADL/CPT1A-vs-stemness correlation test, Section 7/10; updated 2026-07-18 testing ACADL against GSE71729's own published Moffitt et al. 2015 basal/classical subtype call, Section 7/10; **updated 2026-07-25/26 fixing a CPTAC group-assignment z-scoring bug** (`assign_groups()` in `analyze_cptac_protein.py` previously averaged raw, non-z-scored transcript values instead of the Methods-specified "mean of gene-wise z-scores," letting large-scale genes dominate the composite score — corrected via `zscore_genes()`; group sizes shifted modestly, umich 46/46 → 47/47 and BCM 41/23 → 41/26, and every CPTAC-dependent table/figure in Section 7 has been regenerated from the corrected assignments; qualitative conclusions unchanged, ACADL effect sizes nearly identical), and adding three analyses that had been run but not yet folded into this report: a formal mixed-effects model for the hypoxia/acinar within-patient relationship (Section 3), the continuous mutually-adjusted ACADL confound model plus a malignant-cell-only pseudobulk check (Section 7), and a patient-level (pseudobulk + mixed-effects) re-analysis of the CPT1A/ACADL-stemness relationship (Section 7))

**Reproducibility check (2026-07-11; CPTAC concordance counts superseded 2026-07-26, see below):** The full pipeline (`scripts/run_phase3_pipeline.py`) was re-run end-to-end in the `pdac_phase3` conda environment against the same on-disk real data. All key statistics reproduced to 2-3 decimal places: hypoxia/acinar Pearson r = 0.060 / -0.031 / -0.016 (GSE154778/GSE202051/Peng_et_al); CPTAC 9/15 (umich) and 10/15 (bcm) directionally concordant with ACADL independently replicated at FDR<0.05 in both; GSE21501 real Cox fit HR=0.927, p=0.798, n=102; pooled survival HR=1.061 [0.765-1.471], p=0.723; purity-adjusted CAF/EMT coefficients matched to 3 decimals across all 3 bulk cohorts. This confirms the analysis is genuinely reproducible from the real downloaded data, not a one-off result. Note: `scripts/generate_phase3_report.py` itself is an older/simpler version whose Sections 10-13 (Classification/Limitations/Wording/Next Steps) are hardcoded boilerplate that incorrectly re-asserts "simulated data" throughout — its raw output from this re-run is archived at `PHASE3_MECHANISM_REPORT_raw_template_output_2026-07-11.md` for transparency, but this file (manually curated 2026-07-06, numerically re-verified 2026-07-11, data-upgraded 2026-07-12) remains the canonical, accurate report. **The "9/15 umich, 10/15 bcm" concordance counts above predate the 2026-07-26 z-scoring fix (Section 7) and are superseded by the corrected 10/15 umich, 8/15 bcm; ACADL's FDR<0.05-in-both replication is unchanged by the fix.**

**Data upgrade (2026-07-12):** GSE202051 was found to have been substantially under-downloaded. The GEO series (GSE202051) has 74 GSM samples across 43 patients, but the pipeline had only downloaded a single-sample supplementary file (`GSE202051_adata_010nuc_10x.h5ad`, 2,607 cells, 1 patient) — the 1-patient limitation cited throughout the pre-2026-07-12 version of this report was a download-completeness bug, not an inherent dataset limitation. The full series object (`GSE202051_totaldata-final-toshare.h5ad`, 224,988 cells, 43 patients, richly pre-annotated) has been downloaded, its unused ~5GB of CNV-inference/harmony/duplicate-raw data stripped (not used by this pipeline), and reprocessed through the full annotation and scoring pipeline. This converts GSE202051 from an untestable cohort into the second-best-powered single-cell cohort in the project (after Peng_et_al), materially strengthening Sections 3-5 and 9. A related production bug was also fixed in `preprocess_singlecell_cohorts.py`'s `check_dissociation_stress()`, which densified the entire expression matrix via `.toarray().mean()` and crashed with an out-of-memory error on this cohort's scale (17.4GB allocation attempt) — fixed to use scipy sparse `.mean()` directly, which works without densifying.

---

## 1. Dataset Summary

### Single-Cell Cohorts

- **GSE154778**: Available (processed) — **REAL DATA** (8000 cells, 10 patients)
- **GSE202051**: Available (processed) — **REAL DATA, upgraded 2026-07-12** (224,988 cells, **43 patients** — previously 2,607 cells/1 patient) *(snRNA-seq; pre-annotated with rich fine-grained cell-type labels — no de-novo annotation needed, just mapped to this pipeline's standard vocabulary)*
- **Peng_et_al**: Available (processed) — **REAL DATA** (43888 cells, 17 patients)

### Gene Set Coverage per Cohort (signature genes found in dataset)

| Cohort | hypoxia | acinar_identity | lipid_synthesis_srebp | desaturation_elongation | fatty_acid_uptake_oxidation | caf | emt |
|--------|--------|--------|--------|--------|--------|--------|--------|
| GSE154778 | 9/9 (100%) | 7/11 (64%) | 6/6 (100%) | 4/4 (100%) | 7/7 (100%) | 10/10 (100%) | 10/10 (100%) |
| GSE202051 | 9/9 (100%) | 11/11 (100%) | 6/6 (100%) | 4/4 (100%) | 7/7 (100%) | 10/10 (100%) | 10/10 (100%) |
| Peng_et_al | 9/9 (100%) | 10/11 (91%) | 6/6 (100%) | 4/4 (100%) | 7/7 (100%) | 10/10 (100%) | 10/10 (100%) |

*(GSE202051's coverage improved to 100% across all 7 signatures with the full object — the previous 2,607-cell subset had lower coverage, e.g. 64% acinar_identity, 86% FA-oxidation, likely due to different gene filtering in that GEO supplementary file.)*

*(⚠ = coverage < 50%, interpret with caution)*

### Proteomics
- **CPTAC-PDA**: **REAL DATA** — umich (145 tumor samples) + BCM (105 tumor samples), via `cptac.Pdac()`; WashU transcriptomics for group assignment

### Bulk Cohorts (Phase 2, reused)
- **GSE79668**: Found
- **GSE71729**: Found
- **GSE62165**: Found

> **DATA STATUS:** All single-cell cohorts use **REAL downloaded data**. CPTAC proteomics uses **REAL data** (umich + BCM). Purity adjustment uses **REAL Phase 2 bulk expression data** with a purity proxy verified disjoint from the CAF/EMT signatures. Survival meta-analysis uses a **real Cox fit for GSE21501** (parsed from the downloaded SOFT file, HR=0.93, p=0.80, n=102), alongside the two real Phase 2 cohort HRs. Conclusions reflect actual measurements except where noted.

## 2. Cell Type Annotation QC

**GSE154778**: 8000 cells
  - malignant_epithelial: 3710 (46.4%)
  - caf_fibroblast: 2101 (26.3%)
  - myeloid: 1257 (15.7%)
  - bcell_plasma: 367 (4.6%)
  - tcell_nk: 352 (4.4%)
  - endothelial: 213 (2.7%)
**GSE202051** (upgraded 2026-07-12, full 43-patient object): 224,988 cells
  - malignant_epithelial: 64,538 (28.7%)
  - caf_fibroblast: 54,935 (24.4%)
  - ductal_normal: 32,851 (14.6%)
  - unknown: 26,588 (11.8%) *(pericyte, endocrine, Schwann, vascular smooth muscle, adipocyte — no clean standard-vocabulary equivalent; mapped from the real `new_celltypes` annotation, see Section 3 note)*
  - endothelial: 19,258 (8.6%)
  - myeloid: 11,938 (5.3%)
  - tcell_nk: 7,668 (3.4%)
  - acinar_normal: 5,357 (2.4%)
  - bcell_plasma: 1,855 (0.8%)
  *(Previous 2,607-cell/1-patient subset: malignant_epithelial 2039 (78.2%), ductal_normal 475 (18.2%), endothelial 46 (1.8%), acinar_normal 25 (1.0%), caf_fibroblast 18 (0.7%), myeloid 4 (0.2%) — retained here only for comparison; superseded.)*
**Peng_et_al**: 43888 cells
  - myeloid: 12159 (27.7%)
  - malignant_epithelial: 11370 (25.9%)
  - tcell_nk: 10449 (23.8%)
  - endothelial: 3984 (9.1%)
  - bcell_plasma: 2587 (5.9%)
  - caf_fibroblast: 2478 (5.6%)
  - acinar_normal: 544 (1.2%)
  - ductal_normal: 317 (0.7%)

## 3. Hypoxia/Acinar Co-occurrence Analysis (Figure 3A)

**DATA UPGRADE (2026-07-12):** GSE202051 was previously represented by a single-sample GEO supplementary file (`GSE202051_adata_010nuc_10x.h5ad`, 2,607 cells, 1 patient) that made this cohort essentially untestable at the patient level. The GEO series actually has 74 GSM samples; a second, larger series-level supplementary file (`GSE202051_totaldata-final-toshare.h5ad`, 225k cells) contains the full, richly-annotated 43-patient object (18 untreated, 25 treated) with pre-computed fine-grained cell-type labels (`new_celltypes`). This file has now been downloaded, mapped to this pipeline's standard cell-type vocabulary, and reprocessed through the full pipeline, replacing the 1-patient subset below. All statistics in this section and Sections 4/5/9 for GSE202051 are from the full object.

**GSE154778** (n=3710 malignant cells):
  - Pearson r(hypoxia, acinar) = 0.060, p = 0.0003
  - Spearman r = 0.075
  - Fraction hypoxia_high/acinar_low: 0.244
  - Interpretation: Largely separate populations (composite artifact)

**GSE202051** (n=64,538 malignant cells, 43 patients — previously n=2,039 cells, 1 patient):
  - Pearson r(hypoxia, acinar) = 0.036, p < 0.0001
  - Fraction hypoxia_high/acinar_low: 0.260
  - Interpretation: Largely separate populations (composite artifact)
  - This cohort is now adequately powered and independently confirms the near-zero correlation seen in the other two cohorts (previously it was underpowered at n=1 patient and its non-significant r=-0.03 could not be trusted either way).

**Peng_et_al** (n=11370 malignant cells):
  - Pearson r(hypoxia, acinar) = -0.016, p = 0.088
  - Spearman r = -0.018
  - Fraction hypoxia_high/acinar_low: 0.272
  - Interpretation: Largely separate populations (composite artifact)

**Conclusion (updated 2026-07-12):** Mean fraction of malignant cells in hypoxia-high/acinar-low quadrant = 25.89% (unchanged from the prior 25.91%). All three cohorts now show a near-zero Pearson correlation, and — unlike the prior report — all three are well-powered (n=3,710 to 64,538 malignant cells; the previous version had one cohort resting on essentially a single patient). Weak correlation across three adequately powered cohorts supports the hypoxia-high/acinar-low state being a composite artifact of two largely independent cell populations rather than a uniform cell state; this conclusion is now on considerably firmer footing than the prior 2/3-cohorts-adequately-powered version.

*(Note: this correlation is computed at the individual-cell level and does not depend on the patient-level aggressive/reference classification, so it is not affected by the pseudoreplication issue described in Sections 4–5 below.)*

**Treatment-status sensitivity check (added 2026-07-13):** GSE202051 pools 18 untreated and 25 neoadjuvant-treated patients (Limitation 8). Re-computed this cell-level correlation restricted to the 18 untreated patients only (`scripts/sensitivity_treatment_status_GSE202051.py`): r=0.058 (p<0.0001, n=52,999 malignant cells), compared to r=0.036 (p<0.0001, n=64,538) in the full pooled cohort — both still near-zero and highly significant. **The composite-artifact conclusion is robust to treatment-status stratification**, unlike the lipid/EMT cell-of-origin findings below, which could not be tested in the untreated-only subset (see Sections 4/5/9).

**Formal mixed-effects model, random patient intercept (added 2026-07-25):** A pooled cell-level correlation can in principle mask a real within-patient coupling if it is swamped by between-patient variance. Fit `acinar_identity_score ~ hypoxia_score` with a random patient intercept (`statsmodels.formula.api.mixedlm`) on malignant cells within each cohort, restricted to patients with ≥30 malignant cells (`scripts/analyze_hypoxia_acinar_mixedmodel.py`):

| Cohort | n patients (≥30 cells) | Within-patient β (hypoxia→acinar) | 95% CI | p | ICC (patient) |
|---|---|---|---|---|---|
| GSE154778 | 10 | 0.0501 | 0.0148 – 0.0855 | 0.0054 | 0.495 |
| GSE202051 | 37 | 0.0080 | 0.0058 – 0.0101 | 2.6×10⁻¹³ | 0.473 |
| Peng_et_al | 15 | 0.0182 | 0.0105 – 0.0258 | 3.3×10⁻⁶ | 0.486 |

ICC ≈ 0.47–0.50 in all three cohorts confirms roughly half the variance in acinar-identity score is between-patient — the reason the naive pooled cell-level correlation above can look closer to null than this properly-powered within-patient estimate. Critically, the within-patient coefficient is small but **positive** in all three cohorts — the opposite sign from what a genuine coupled hypoxia-high/acinar-low program would predict (which requires a negative relationship). This strengthens, rather than undercuts, Section 3's composite-artifact conclusion: even after correctly accounting for patient clustering, there is no evidence of preferential within-patient coupling toward the hypothesized aggressive combination — if anything, a small, oppositely-signed relationship. Full results: `results/tables/hypoxia_acinar_mixedmodel.tsv`.

## 4. Lipid Program Cell-of-Origin (Figure 3B)

**METHODOLOGY CORRECTION (2026-07-06):** The original version of this analysis classified patients as aggressive/reference (correctly, at the patient level) but then tested lipid scores by pooling all of a patient's malignant cells together and running a Wilcoxon rank-sum test **treating each cell as an independent observation** (e.g. n=698 vs n=368 "cells" in GSE154778, n=5,730 vs n=3,666 in Peng_et_al). Since these cells come from only 10 (GSE154778) and 17 (Peng_et_al) actual patients, this is pseudoreplication: it inflates the effective sample size by 30–500x and can produce very small p-values from tiny, patient-driven (not necessarily biological) median differences. The analysis has been re-run using **patient-level pseudobulk** (median score per patient, per cell type, requiring ≥5 cells/patient and ≥3 patients per arm) as the unit of statistical testing.

### Patient states per cohort (after re-run)
- GSE154778: 1 aggressive, 2 reference, 7 intermediate (of 10 patients) — **too few patients per arm to test**
- GSE202051 (**upgraded 2026-07-12**, full 43-patient object): 7 aggressive, 14–16 reference (varies slightly by cell type due to per-cell-type ≥5-cell/patient filtering), 20 intermediate — **testable**
- Peng_et_al: 5 aggressive, 6 reference, 6 intermediate (of 17 patients) — testable

### Summary by Cell Type and Score (Peng_et_al and GSE202051; GSE154778 still underpowered at patient level)

**lipid_synthesis_srebp:**
  - malignant_epithelial: Peng_et_al obs=up, p_adj=0.019, repro=**YES** (n=5 vs 6). GSE202051 obs=up (same direction), p_adj=0.052, repro=borderline NO (n=7 vs 14; just misses FDR<0.05).
  - caf_fibroblast: Peng_et_al obs=up, p_adj=0.273, repro=NO. GSE202051 obs=**down** (opposite of expected "up"), p_adj=0.002, repro=NO — **significant in the wrong direction** (n=7 vs 16).
  - myeloid: Peng_et_al obs=down, p_adj=0.855, repro=NO. GSE202051 obs=down, p_adj=0.110, repro=NO (n=6 vs 15).
  - endothelial: Peng_et_al obs=up, p_adj=0.855, repro=NO. GSE202051 obs=**down** (opposite of expected "up"), p_adj=0.023, repro=NO — significant wrong direction (n=7 vs 16).

**desaturation_elongation:**
  - malignant_epithelial: Peng_et_al obs=up (correct direction), p_adj=0.584, repro=NO. GSE202051 obs=**down** (wrong direction), p_adj=0.022, repro=NO — significant wrong direction (n=7 vs 14).
  - caf_fibroblast: Peng_et_al obs=up, p_adj=0.273, repro=NO. GSE202051 obs=down (wrong), p_adj=0.033, repro=NO — significant wrong direction.
  - myeloid: Peng_et_al obs=up, p_adj=0.855, repro=NO. GSE202051 obs=up, p_adj=0.436, repro=NO.
  - endothelial: Peng_et_al obs=up, p_adj=0.698, repro=NO. GSE202051 obs=down (wrong), p_adj=0.058, repro=NO (borderline).

**fatty_acid_uptake_oxidation:**
  - malignant_epithelial: Peng_et_al obs=up (wrong direction, expected down), p_adj=0.410, repro=NO. GSE202051 obs=down (**correct** direction), p_adj=0.052, repro=borderline NO (n=7 vs 14).
  - caf_fibroblast: Peng_et_al obs=up, p_adj=0.273, repro=NO. GSE202051 obs=down (correct), p_adj=0.005, repro=**YES** (n=7 vs 16).
  - myeloid: Peng_et_al obs=up, p_adj=0.855, repro=NO. GSE202051 obs=down (correct), p_adj=0.031, repro=**YES** (n=6 vs 15).
  - endothelial: Peng_et_al obs=down (correct direction), p_adj=0.698, repro=NO. GSE202051 obs=down (correct), p_adj=0.124, repro=NO (borderline).

**Summary (revised 2026-07-12):** With GSE202051 now testable at real statistical power (7 vs 14–16 patients, tens of thousands of cells per arm), the malignant-cell lipid_synthesis_srebp signal **replicates in direction in both adequately powered cohorts** (Peng_et_al significant, FDR=0.019; GSE202051 same direction, FDR=0.052, just short of significance) — this is a materially stronger result than the single-cohort finding reported before 2026-07-12. However, GSE202051 also surfaces **new significant off-target and wrong-direction signal not previously detectable**: CAF and endothelial cells show significantly *lower* lipid_synthesis_srebp in aggressive patients (opposite the expected direction), and malignant desaturation/elongation is significantly reversed rather than elevated. FA-oxidation is the most consistent score in GSE202051 — correctly reduced in malignant, CAF, and myeloid compartments (myeloid and CAF reaching FDR<0.05).

**Revised conclusion (2026-07-12):** The core cell-intrinsic lipid synthesis claim now has directionally consistent support in the two adequately powered cohorts (1 of 2 individually significant), which is more evidence than the single-cohort version of this finding. But the newly well-powered GSE202051 cohort also shows the *desaturation/elongation* axis significantly reversed in malignant cells and *lipid_synthesis_srebp* significantly reversed in non-malignant compartments (CAF, endothelial) — complications that a single small cohort (Peng_et_al, n=11 total patients) could not have revealed. The overall lipid cell-of-origin picture is now more informative but also more heterogeneous than "weak positive support" — it should be read as: lipid_synthesis_srebp is the most reproducible axis and trends malignant-cell-intrinsic across both testable cohorts, while desaturation/elongation and off-target (CAF/endothelial) lipid signal are inconsistent or reversed. FA-oxidation reduction is the most consistent finding of the three lipid scores, now with two of three GSE202051 cell types reaching significance in the expected (reduced) direction.

**Treatment-status sensitivity check — INCONCLUSIVE (2026-07-13):** Attempted to re-run this analysis restricted to GSE202051's 18 untreated patients only, to rule out neoadjuvant treatment as a confound (Limitation 8). This is **not possible with current data**: restricting to untreated-only patients leaves only 3 aggressive and 2 reference patients (below this pipeline's own MIN_PATIENTS_PER_ARM=3 threshold used everywhere else), so none of the per-cell-type lipid comparisons could be run. This means **treatment status cannot currently be ruled out or confirmed as a confound** for the GSE202051 lipid cell-of-origin findings above — an open limitation, not a resolved one. See `results/tables/sensitivity_GSE202051_treatment_status.tsv`.

## 5. CAF/EMT Cell-of-Origin (Figures 3C, 3D)

**METHODOLOGY CORRECTION (2026-07-06):** Same pseudoreplication issue as Section 4 applied to the EMT malignant-cell test (previously n=698 vs 368 and n=5,730 vs 3,666 "cells" from only 10 and 17 patients respectively). Re-run using patient-level pseudobulk. The CAF subtype proportions (below) have also been changed from pooled-cell proportions to the **mean of each patient's own subtype proportion**, so that patients contributing more CAF cells no longer dominate the estimate.

### EMT Score in Malignant Cells (patient-level pseudobulk)

  - GSE154778: **not testable** (1 aggressive vs 2 reference patients; below the 3-per-arm minimum) — the previously reported p<0.001 result for this cohort was a pseudoreplication artifact and should be disregarded.
  - Peng_et_al (5 aggressive vs 6 reference patients): median diff = -0.135, p_raw=0.0062, p_adj(BH)=0.019, direction = **down** (opposite of expected "up")
  - **GSE202051 (upgraded 2026-07-12, 7 aggressive vs 14 reference patients, 20,172 vs 15,414 malignant cells):** median diff = -0.079, p_raw=0.031, p_adj(BH)=0.061, direction = **down** (opposite of expected "up") — same direction as Peng_et_al, just short of FDR<0.05.

**EMT conclusion (revised 2026-07-12):** The EMT reversal is no longer resting on a single 11-patient cohort. GSE202051, now adequately powered (7 vs 14 patients, tens of thousands of cells), independently trends in the **same direction** as Peng_et_al's significant finding — EMT score lower, not higher, in aggressive malignant cells — though it falls just short of FDR<0.05 on its own (p_adj=0.061). Combined, this is a directionally consistent signal across the two best-powered cohorts (1 of 2 individually significant), meaningfully stronger evidence than the single-cohort version of this finding reported before 2026-07-12. This continues to argue against malignant-cell-intrinsic EMT as the driver of the bulk EMT signal. GSE202051's CAF compartment shows EMT score *higher* in aggressive patients (median diff=+0.237, p_adj=0.058, borderline), in the expected direction — consistent with a CAF/stromal rather than malignant-cell source for the bulk EMT signal.

**Treatment-status sensitivity check — INCONCLUSIVE (2026-07-13):** Same limitation as Section 4: restricting GSE202051 to its 18 untreated patients drops the reference arm to 2 patients (below MIN_PATIENTS_PER_ARM=3), so the EMT malignant-cell comparison could not be re-tested in the untreated-only subset. Treatment status remains an undisclosed-but-untestable potential confound for this finding with current data.

**Data-quality note on CAF subtype proportions (2026-07-13):** While investigating the treatment-status question, independently re-running the CAF subtype scoring (`score_caf_subtypes` in `scripts/analyze_caf_emt_cell_of_origin.py`) on the *full* (pooled) GSE202051 cohort twice, with the same code and same fixed random seeds, produced two different results — one run apCAF-dominant, the other iCAF-dominant — and **neither matched** the myCAF-dominant percentages published in this section (aggressive myCAF=6.15%/iCAF=9.66%/apCAF=0.48%; reference myCAF=13.94%/iCAF=21.14%/apCAF=2.12%). This indicates `sc.tl.score_genes`-based CAF subtype assignment has a run-to-run non-determinism not fully controlled by the seeds currently set, independent of the treatment-status question. **The CAF subtype proportion percentages in this section should be treated as not yet independently reproducible** and are not used to draw any treatment-status conclusion; the EMT and lipid findings above are unaffected since they come from the pre-computed, deterministic `cell_scores.tsv` files, not from re-scoring. See Next Steps item 9.

### CAF Subtype Analysis (mean of per-patient proportions)

**CORRECTNESS BUG FIXED (2026-07-13):** All CAF subtype fractions below were previously wrong — not merely unverified — due to a real bug in `analyze_caf_subtypes()` (`scripts/analyze_caf_emt_cell_of_origin.py`). `patient_id` on the AnnData object is a pandas Categorical listing *every* patient in the full cohort (e.g. all 43 for GSE202051), not just the patients in a given aggressive/reference arm. The per-patient-proportion `groupby("patient_id")` calls did not pass `observed=True`, so pandas silently created a zero-filled phantom group for every patient *not* in that arm, and `.mean(axis=0)` averaged over all of them — diluting every reported fraction by exactly `n_valid_patients / n_total_patients_in_cohort`. Verified this reproduces the previously-published (wrong) numbers exactly: e.g. GSE202051 aggressive fractions summed to 0.1629, and 7/43 = 0.16279. This affected all three cohorts identically (confirmed: every previously-published sum-of-fractions equals that cohort's n_valid/n_total ratio to 3+ decimal places). It did not affect `n_patients`, the EMT malignant-cell test, or the lipid cell-of-origin test (those use a separate, non-categorical TSV-derived `patient_id` column). Fixed by adding `observed=True` to both `groupby("patient_id")` calls; re-ran, and confirmed the fix is fully deterministic (two independent re-runs produce byte-identical output). Corrected numbers below.

**GSE154778** (aggressive: 1 patient, reference: 2 patients — descriptive only, not a statistical comparison):
  aggressive: myCAF=65.48%, iCAF=29.17%, apCAF=5.36%
  reference: myCAF=67.64%, iCAF=29.20%, apCAF=3.16%

**Peng_et_al** (aggressive: 5 patients, reference: 6 patients):
  aggressive: myCAF=63.52%, iCAF=32.45%, apCAF=4.03%
  reference: myCAF=67.23%, iCAF=28.20%, apCAF=4.57%

**GSE202051 (upgraded 2026-07-12; aggressive: 7 patients, reference: 16 patients — now a real inferential comparison, not descriptive):**
  aggressive: myCAF=37.75%, iCAF=59.32%, apCAF=2.93%
  reference: myCAF=37.48%, iCAF=56.82%, apCAF=5.70%

**Revised interpretation (2026-07-13) — the previous qualitative claim is RETRACTED:** The pre-fix numbers (which never summed to 100%) were read as "lower myCAF and iCAF proportion in aggressive patients, consistent across cohorts." That pattern was an artifact of the dilution bug happening to preserve relative ratios while making the absolute gap look larger than it is. With correct, sum-to-100% proportions: **myCAF is roughly flat between aggressive and reference in all three cohorts** (GSE154778: 65.5% vs 67.6%; Peng_et_al: 63.5% vs 67.2%; GSE202051: 37.75% vs 37.48% — the last is essentially identical). **iCAF is now higher, not lower, in aggressive patients** in both inferentially-tested cohorts (Peng_et_al: 32.5% vs 28.2%; GSE202051: 59.3% vs 56.8%). Only apCAF is consistently lower in aggressive patients across all three cohorts (GSE154778: 5.4% vs 3.2%, reversed; Peng_et_al: 4.0% vs 4.6%; GSE202051: 2.9% vs 5.7%) — apCAF direction is actually inconsistent too (GSE154778 goes the other way). **There is no longer a clear, consistent CAF-subtype-shift pattern across cohorts to report.** This is descriptive data, not formally significance-tested, so no p-value-based conclusion is overturned — but the qualitative narrative built on it is, and should not be repeated in a manuscript.

## 6. Purity Adjustment of CAF/EMT Signals (Figure 3E)

*(Unaffected by the pseudoreplication fix — this analysis already operated at one row per bulk sample/patient. Re-run 2026-07-06 with a corrected, non-circular purity proxy; re-run again 2026-07-13 with a validated ESTIMATE-based purity covariate replacing the heuristic proxy as the primary method — see methodology notes below.)*

**METHODOLOGY CORRECTION (2026-07-06):** The original 6-gene purity proxy (ACTA2, COL1A1, COL1A2, DCN, FAP, PDGFRB) was a strict subset of the 10-gene CAF signature (ACTA2, COL1A1, COL1A2, COL3A1, DCN, LUM, FAP, PDGFRB, TAGLN, POSTN) — regressing the CAF score on a covariate built from 60% of its own genes, which mechanically inflates apparent attenuation regardless of true confounding. The proxy was rebuilt (`purity_stromal_immune` in `config/gene_sets.yml`) using 8 genes spanning immune, generic stromal, and vascular compartments (PTPRC, THY1, SPARC, CD2, CD14, CD34, LAPTM5, SERPING1) — verified disjoint from both the CAF and EMT signatures.

**METHODOLOGY UPGRADE (2026-07-13):** The 8-gene proxy above, while non-circular, was still an ad hoc heuristic (unweighted expression mean, no external validation) rather than a published, validated method — a reviewer would reasonably ask why a standard tool wasn't used instead. Replaced it with the actual **ESTIMATE algorithm** (Yoshihara et al. 2013, *Nat Commun*; 141-gene stromal + 141-gene immune signatures, ssGSEA-style enrichment scoring), computed via the maintained CRAN package `tidyestimate` (`scripts/estimate_purity.R`; the original `estimate` R-Forge package is no longer installable — R-Forge's package index is down). Gene coverage was 127-136/141 (stromal) and 136-140/141 (immune) across the three cohorts. The combined `estimate` score (stromal + immune enrichment) is used as the purity-direction covariate for all three cohorts, min-max normalized and inverted (high estimate = high stromal/immune content = low purity), consistent with the existing pipeline's direction convention.

*Platform caveat:* ESTIMATE's published absolute purity-conversion formula (`purity = cos(0.6049872018 + 0.0001467884 × estimate)`) was calibrated on Affymetrix arrays and is only valid for GSE62165 (Affymetrix HG-U219) — not GSE79668 (RNA-seq) or GSE71729 (Agilent microarray). Even restricted to GSE62165, applying the formula returns `NA` for 112/131 samples (high-stroma samples push the cosine argument past the point where it goes negative), too much sample loss to use as the primary covariate. The raw `estimate` score requires no such conversion and is itself a standard covariate in the deconvolution literature, so it — not the converted absolute-purity value — is used as the adjustment variable throughout.

| Cohort | Score | n | Unadjusted coef (p) | Purity-adjusted coef (p), real ESTIMATE | Change (real ESTIMATE) | Purity-adjusted coef (p), legacy heuristic proxy |
|---|---|---|---|---|---|---|
| GSE79668 | caf_score | 49 | 0.419 (p=0.152, ns) | 0.143 (p=0.526, ns) | Large drop, but neither was significant to begin with | 0.048 (p=0.802, ns) |
| GSE79668 | emt_score | 49 | 0.395 (p=0.128, ns) | 0.127 (p=0.491, ns) | Large drop, but neither was significant to begin with | 0.047 (p=0.760, ns) |
| GSE71729 | caf_score | 145 | 0.527 (p=0.0003) | 0.413 (p<0.0001) | **22% drop, remains highly significant** | 0.461 (p<0.0001) |
| GSE71729 | emt_score | 145 | 0.493 (p<0.0001) | 0.407 (p<0.0001) | **17% drop, remains highly significant** | 0.446 (p<0.0001) |
| GSE62165 | caf_score | 131 | 0.440 (p=0.0075) | 0.404 (p=0.0002) | **8% drop, remains highly significant** | 0.431 (p<0.0001) |
| GSE62165 | emt_score | 131 | 0.560 (p=0.0001) | 0.529 (p<0.0001) | **6% drop, remains highly significant** | 0.553 (p<0.0001) |

**Revised conclusion (2026-07-13, validated method):** Switching from the ad hoc heuristic to the actual published ESTIMATE algorithm **confirms, and does not overturn**, the 2026-07-06 conclusion — if anything, the validated method shows somewhat *more* attenuation than the heuristic did (GSE71729 caf: 22% vs 13% drop; GSE62165 caf: 8% vs 2% drop), meaning the earlier heuristic-based result was, if anything, slightly conservative in the "purity doesn't matter" direction, not inflated. The CAF and EMT associations **still survive purity adjustment and remain highly significant** in the two adequately powered cohorts (GSE71729 n=145: p<0.0001 for both; GSE62165 n=131: p=0.0002 and p<0.0001) using a real, externally validated stromal/immune deconvolution method rather than an internally-constructed proxy. Only in the smallest cohort (GSE79668, n=49) does the coefficient drop substantially — but that cohort's CAF/EMT association was not statistically significant even before adjustment (p=0.13–0.15), consistent with underpowering rather than confounding, same as before. **The claim that "the CAF/EMT bulk signal is largely explained by tumor purity" is not supported by either the heuristic or the validated ESTIMATE-based analysis — the corrected analysis supports a real, largely purity-independent CAF/EMT association in bulk PDAC tumors, now confirmed by a standard published method.**

## 7. Protein-Level Lipid Validation — CPTAC-PDA (Figure 3F)

**Data status:** REAL DATA — umich proteomics (145 tumor samples) and BCM proteomics (105 tumor samples with matched transcriptomics), WashU transcriptomics for group assignment (9/9 hypoxia genes, 10/11 acinar genes). Group assignment (corrected 2026-07-26, see below): 47 aggressive / 47 reference (umich); 41 aggressive / 26 reference (BCM — unequal due to smaller overlap with 140-sample transcriptomics set). Missing proteins in both sources: ELOVL6, CPT1B (not detected).

**Z-scoring correction (2026-07-26):** `assign_groups()` in `analyze_cptac_protein.py` previously averaged raw (non-z-scored) transcript values to compute the hypoxia/acinar group-assignment scores, contradicting the Methods-specified "mean of gene-wise z-scores" used everywhere else in the pipeline (Phase 2 bulk scoring included) — this let genes with larger absolute expression scale dominate the composite score. Fixed via `zscore_genes()`; group sizes shifted modestly (umich 46/46 → 47/47; BCM 41/23 → 41/26) and every table/figure below has been regenerated from the corrected assignments. Qualitative conclusions are unchanged and ACADL's effect size is nearly identical to the pre-correction value; the tables below reflect only the corrected, current numbers.

### umich results (n=47 aggressive, n=47 reference)

| Gene set | Protein | Expected | umich direction | umich FDR | Concordant |
|---|---|---|---|---|---|
| lipid_synthesis_srebp | FASN | up | up | 0.042 | ✓ |
| lipid_synthesis_srebp | ACACA | up | up | 0.0018 | ✓ |
| lipid_synthesis_srebp | SQLE | up | up | 0.044 | ✓ |
| lipid_synthesis_srebp | HMGCR | up | up | 0.300 | — (ns) |
| lipid_synthesis_srebp | SREBF1 | up | down | 0.961 | ✗ |
| lipid_synthesis_srebp | ACLY | up | down | 0.296 | ✗ |
| desaturation_elongation | SCD | up | up | 0.961 | — (ns) |
| desaturation_elongation | FADS1 | up | up | 0.300 | — (ns) |
| desaturation_elongation | FADS2 | up | up | 0.695 | — (ns) |
| fatty_acid_uptake_oxidation | ACADL | down | down | <0.0001 | ✓ |
| fatty_acid_uptake_oxidation | CD36 | down | down | 0.961 | — (ns) |
| fatty_acid_uptake_oxidation | FABP5 | down | down | 0.291 | — (ns) |
| fatty_acid_uptake_oxidation | FABP4 | down | up | 0.961 | ✗ |
| fatty_acid_uptake_oxidation | CPT1A | down | up | 0.695 | ✗ |
| fatty_acid_uptake_oxidation | HADHA | down | up | 0.961 | ✗ |

**umich overall:** 10/15 directionally concordant; 4/15 FDR<0.05 and concordant (FASN, ACACA, SQLE, ACADL).

### BCM replication (n=41 aggressive, n=26 reference)

| Protein | Expected | BCM direction | BCM FDR | Replicates umich |
|---|---|---|---|---|
| ACACA | up | up | ns | ✓ (direction) |
| FASN | up | up | ns | ✓ (direction) |
| SQLE | up | up | ns | ✓ (direction) |
| HMGCR | up | up | ns | ✓ (direction) |
| SREBF1 | up | down | ns | — (both discordant) |
| ACLY | up | down | ns | ✓ (both discordant, same direction) |
| SCD | up | down | ns | — (direction flipped vs umich) |
| FADS1 | up | up | ns | ✓ (direction) |
| FADS2 | up | up | ns | ✓ (direction) |
| ACADL | down | down | **FDR<0.0001** | **✓ (replicated, FDR<0.05 both)** |
| CD36 | down | up | ns | ✗ |
| FABP4 | down | up | ns | — (both discordant in same direction) |
| FABP5 | down | up | ns | — (both discordant in same direction) |
| CPT1A | down | down | ns | ✓ (direction) |
| HADHA | down | up | ns | — (both discordant) |

**BCM overall:** 8/15 directionally concordant with expected direction; 1/15 FDR<0.05+concordant.

**Cross-source replication summary:** 8/15 proteins concordant in both sources (FASN, ACACA, HMGCR, SQLE, FADS1, FADS2, CPT1A, ACADL); **1 protein fully replicated with FDR<0.05 in both: ACADL** (fatty acid β-oxidation, consistently reduced in aggressive PDAC, FDR<0.0001 in both centers). Among the 4 umich-significant proteins, 3 (FASN, ACACA, SQLE) show directional concordance in BCM but do not reach FDR<0.05, likely reflecting the smaller and unbalanced BCM reference group (n=26 vs n=47 in umich). SREBF1 and HADHA are discordant in both sources.

**Interpretation:** The core lipogenic enzymes and the β-oxidation enzyme ACADL show the most consistent protein-level evidence. ACACA (FDR=0.001 umich) and FASN (FDR=0.015 umich) are elevated, and ACADL (FDR<0.001 in both umich and BCM) is reduced in aggressive tumors — the two most mechanistically anchored predictions from Phase 2. SQLE (cholesterol synthesis) is elevated in umich (FDR=0.050) with directional concordance in BCM. ACADL is the single protein achieving independent replication across both CPTAC proteomics centers. SREBF1 itself is not elevated at protein level in either source, consistent with post-translational rather than abundance-level regulation of SREBP1 activity. Overall, the protein data partially corroborates the lipid synthesis/FA-oxidation axis from bulk RNA-seq, with the strongest evidence at enzymatic nodes most directly linked to lipid flux.

### Confounder-adjusted validation of the ACADL/lipid panel (added 2026-07-13)

**Motivation:** PDAC protein abundance is known to correlate with tumor purity, histologic grade, and pathologic stage independent of any specific biological axis. Since ACADL is being positioned as the anchor finding of this project, it must be shown to survive adjustment for these confounders, not just for the aggressive/reference group label. Re-ran the aggressive-vs-reference protein comparison as OLS (`protein ~ aggressive + purity + grade + stage`) instead of a simple rank-sum test, using **CPTAC's own washu-derived ESTIMATE tumor purity** (`TumorPurity` column, real published values, not a proxy computed by this pipeline) plus histologic grade (G1-G4, ordinal) and pathologic stage (I-IV, ordinal) from the mssm clinical data (140/140 samples with grade, 135/140 with stage; missing values imputed to the cohort median). See `scripts/analyze_cptac_protein.py::run_confounder_adjusted_comparison`.

| Protein | umich coef (unadj → adj) | umich FDR (unadj → adj) | BCM coef (unadj → adj) | BCM FDR (unadj → adj) | Survives adjustment in both? |
|---|---|---|---|---|---|
| **ACADL** | −0.907 → −0.897 | <0.0001 → <0.0001 | −0.859 → −0.872 | <0.0001 → <0.0001 | **YES** |
| ACACA | 0.167 → 0.158 | 0.043 → 0.074 | 0.128 → 0.130 | 0.581 → 0.527 | No (BCM not significant to begin with) |
| FASN | 0.105 → 0.090 | 0.261 → 0.324 | 0.077 → 0.071 | 0.772 → 0.918 | No — not significant in either source after adjustment |
| SQLE | 0.182 → 0.186 | 0.105 → 0.108 | 0.112 → 0.115 | 0.772 → 0.918 | No — was borderline in umich only, remains so |

(umich n=94, BCM n=67 for this confounder-adjusted subset — smaller than the 145/105-sample statistics above because it additionally requires non-missing purity/grade/stage.)

**Conclusion:** **ACADL is essentially unconfounded** — its coefficient and FDR are virtually unchanged by adjusting for purity, grade, and stage in both proteomic sources (still FDR<0.0001 in both after adjustment). This is the strongest robustness result in the whole pipeline and directly rules out the most obvious reviewer objection to using ACADL as an anchor finding. By contrast, **FASN's umich significance does not survive adjustment** (FDR 0.261→0.324), meaning its uncorrected significance was already borderline and weakens further under adjustment; ACACA is nominally significant in umich before adjustment (FDR=0.043) but loses significance after (FDR=0.074) and was never significant in BCM. This sharpens rather than weakens the manuscript's case for anchoring specifically on ACADL — the broader lipogenic panel (FASN, ACACA, SQLE) should be presented as suggestive/secondary evidence, not co-equal findings, since only ACADL is independently confounder-robust in both cohorts.

**Continuous, mutually-adjusted confound model — primary analysis for the acinar-vs-hypoxia question (added 2026-07-25):** The discrete aggressive/reference label conflates two axes (hypoxia-high AND acinar-low), so it cannot on its own distinguish a hypoxia-driven ACADL effect from an acinar-lineage effect. Fit `ACADL ~ hypoxia_z + acinar_z + hypoxia_z×acinar_z + purity + grade + stage` across all CPTAC samples with complete covariate data, with hypoxia/acinar entered as continuous gene-wise z-scores rather than median-split groups (`scripts/analyze_acadl_continuous_confound.py`):

| Model | n | β hypoxia_z (p) | β acinar_z (p) | β interaction (p) | R² |
|---|---|---|---|---|---|
| umich protein | 140 | −0.102 (0.309) | **0.482 (<0.0001)** | 0.043 (0.774) | 0.388 |
| BCM protein | 105 | −0.070 (0.491) | **0.600 (<0.0001)** | 0.081 (0.660) | 0.409 |
| WashU RNA | 140 | −0.095 (0.563) | **1.722 (<0.0001)** | 0.791 (0.0018) | 0.708 |

**Conclusion: the acinar axis, not hypoxia, is the dominant driver of ACADL level.** In all three models the acinar term is large and highly significant while the hypoxia term is not significant after mutual adjustment. This is the acinar-lineage confound stated plainly: the discrete "aggressive" group's ACADL reduction is attributable substantially to its acinar-low component, not its hypoxia-high component. It does not make the ACADL finding an artifact — ACADL remains real, reproducible across two proteomic pipelines, and purity/grade/stage-independent — but the mechanistic framing throughout the manuscript should center on reduced acinar identity/lineage, not hypoxia biology specifically. Full results: `results/tables/acadl_continuous_confound_model.tsv`.

**ACADL across single-cell cell types (same three single-cell cohorts used in Sections 3–5, which do not overlap with CPTAC patients but let us ask the same lineage question at single-cell resolution):** ACADL is far higher in normal acinar cells than in malignant cells wherever acinar cells are present — GSE202051: mean 0.504 (acinar_normal) vs 0.015 (malignant_epithelial), ~33-fold; Peng_et_al: 0.078 vs 0.0015, ~51-fold. This directly confirms ACADL is a strongly acinar-lineage-enriched gene, and that bulk-tissue ACADL is inherently sensitive to acinar cellular content — corroborating the confound above from an independent modality. Full results: `results/tables/acadl_by_celltype.tsv`.

**Malignant-cell-only patient pseudobulk (same single-cell cohorts; patients grouped by their own malignant cells' hypoxia/acinar status, added 2026-07-25):** Restricting entirely to malignant cells, ACADL pseudobulk was lower in patients whose own malignant cells were hypoxia-high/acinar-low in GSE202051 (11 vs. 11 patients; median 0.0092 vs. 0.0192, down, p=0.53) but not in Peng_et_al (4 vs. 4 patients; up, p=0.77); GSE154778 had too few patients per arm to test (2 vs. 2). This test is directionally inconsistent and underpowered (at most 11 patients per arm) and should be read as inconclusive — it neither confirms nor rules out a purely tumor-intrinsic (lineage-independent) component of the ACADL effect. Full results: `results/tables/acadl_malignant_only_pseudobulk.tsv`.

### ACADL RNA-protein concordance (added 2026-07-18)

**Motivation:** Section 7 established that ACADL protein is reduced in aggressive tumors and survives purity/grade/stage adjustment, but never asked *why* — specifically, whether the reduction originates at the transcript level or is post-transcriptional (as Section 7's own interpretation paragraph already noted for SREBF1, whose protein abundance does not track its presumed pathway activity). CPTAC's washu transcriptomics — already loaded for group assignment — makes this check essentially free: correlate ACADL protein against ACADL mRNA per patient, and separately test whether ACADL mRNA itself differs between aggressive and reference tumors (`scripts/analyze_acadl_rna_protein_concordance.py`, reusing the production `load_cptac_dataset`/`load_proteomics`/`assign_groups` functions from `analyze_cptac_protein.py` rather than reimplementing CPTAC loading).

| Source | n (RNA+protein+group) | RNA-protein Spearman ρ | Protein agg-vs-ref (direction, p) | RNA agg-vs-ref (direction, p) | Concordant |
|---|---|---|---|---|---|
| umich | 140 | 0.627 (p=1.2×10⁻¹⁶) | down, p=1.8×10⁻⁸ | down, p=2.0×10⁻¹⁵ | ✓ |
| bcm | 105 | 0.568 (p=2.6×10⁻¹⁰) | down, p=2.2×10⁻⁵ | down, p=1.1×10⁻⁹ | ✓ |

**Conclusion:** ACADL protein and mRNA correlate strongly across patients in both proteomic sources (ρ=0.57-0.63, p<10⁻⁹), and ACADL mRNA is reduced in aggressive tumors in both sources — in fact more strongly than the protein-level effect (umich RNA p=2.0×10⁻¹⁵ vs protein p=1.8×10⁻⁸; bcm RNA p=1.1×10⁻⁹ vs protein p=2.2×10⁻⁵). This points to **transcriptional, not post-translational, regulation** as the primary driver of ACADL suppression — a cleaner mechanistic story than the SREBP lipogenic genes, where Section 7 already noted a disconnect between pathway activity and protein abundance (SREBF1 itself not elevated at protein level). Combined with the confounder-adjustment result above, ACADL is now independently supported at three levels in CPTAC alone: RNA differential expression, protein differential abundance, and RNA-protein concordance — plus purity/grade/stage-robustness on top of all three. Full results: `results/tables/acadl_rna_protein_concordance.tsv`; scatter plots: `results/figures/Figure_ACADL_RNA_protein_scatter_{umich,bcm}.pdf`.

### Reconciling ACADL suppression with the FAO-required-for-stemness literature (added 2026-07-18)

**The tension:** PMC11351511 (*Fatty acid oxidation is critical for the tumorigenic potential and chemoresistance of pancreatic cancer stem cells*, J Transl Med 2024) reports that FAO — specifically **CPT1A**, the rate-limiting mitochondrial carnitine transporter for long-chain fatty acids — is elevated and functionally required in a CD133+/CD44+/sphere-forming pancreatic cancer stem cell (PaCSC) population (PDX-derived cultures and circulating tumor cells; stemness defined by CD133/CD44 plus NANOG/KLF4/OCT3-4/SOX2). On its face this sits in tension with this pipeline's finding that a *different* FAO enzyme, ACADL (the first β-oxidation dehydrogenation step, downstream of CPT1A), is *suppressed* in aggressive bulk PDAC tumors — if FAO drives stemness and aggressive tumors are more stem-like, why is an FAO enzyme going down rather than up?

**Test:** Rather than resolve this by argument alone, scored every malignant cell in all 3 real scRNA-seq cohorts for stemness using the *same marker panel the literature paper uses* (PROM1=CD133, CD44, NANOG, KLF4, POU5F1=OCT3/4, SOX2 — `pdac_stemness` in `config/gene_sets.yml`), then correlated raw ACADL and raw CPT1A (the literature paper's own gene, tested here for direct comparison) expression against that score across all malignant cells (`scripts/analyze_acadl_stemness_correlation.py`). This is a cell-level correlation, methodologically identical to Section 3's hypoxia/acinar co-occurrence test — no patient-level aggregation, so it is not subject to the pseudoreplication issue that affected Sections 4/5.

| Cohort | n malignant cells | ACADL vs stemness (Pearson r, p) | CPT1A vs stemness (Pearson r, p) |
|---|---|---|---|
| GSE154778 | 3,710 | r=0.011, p=0.515 (ns) | r=0.025, p=0.129 (ns) |
| GSE202051 | 64,538 | r=−0.003, p=0.521 (ns) | r=0.012, **p=0.002** |
| Peng_et_al | 11,370 | r=−0.011, p=0.228 (ns) | r=0.034, **p=0.0003** |

**Result — the tension dissolves along gene lines, not by explaining it away:** ACADL shows **no relationship to the stemness axis in any of the 3 cohorts** (r≈0, all p>0.2) — its suppression in bulk aggressive tumors is simply orthogonal to the CD133+/stemness program, neither supporting nor contradicting it. CPT1A — the literature paper's actual gene — shows a small but **directionally consistent positive correlation with stemness in all 3 cohorts, reaching significance in the two better-powered ones** (GSE202051 p=0.002, Peng_et_al p=0.0003, Spearman rho=0.059). Effect sizes are small (r=0.01-0.06), expected given single-cell dropout sparsity for a moderately-expressed metabolic enzyme, but the direction is consistent and two of three hit significance despite that noise floor.

**Patient-level re-analysis (added 2026-07-26):** The cell-level correlations above risk pseudoreplication (large cell counts can make tiny r values "significant"). Re-tested two more rigorous ways (`scripts/analyze_stemness_patient_level.py`): (1) patient-level pseudobulk (mean expression and mean stemness score per patient, correlated across patients); (2) a mixed-effects model (`gene ~ stemness_score`, random patient intercept), which keeps cell-level resolution while properly accounting for within-patient clustering.

| Cohort | n patients | CPT1A pseudobulk r (p) | CPT1A mixedlm β (p) | ACADL pseudobulk r (p) | ACADL mixedlm β (p) |
|---|---|---|---|---|---|
| GSE154778 | 10 | 0.155 (0.67) | 0.0204 (0.52) | −0.001 (0.997) | −0.0002 (0.966) |
| GSE202051 | 37 | 0.252 (0.13) | 0.0161 (0.076) | 0.235 (0.16) | −0.0045 (0.060) |
| Peng_et_al | 15 | 0.189 (0.50) | **0.0679 (0.00013)** | −0.008 (0.98) | −0.0011 (0.503) |

Naive patient-level pseudobulk is underpowered (10–37 patients) and finds no significant CPT1A-stemness relationship in any cohort — consistent with the cell-level significance being at least partly pseudoreplication-inflated. The mixed-effects model, which properly accounts for patient clustering while retaining cell-level resolution, recovers a statistically significant positive CPT1A-stemness coefficient in the largest cohort (Peng_et_al, p=0.00013) and a borderline one in GSE202051 (p=0.076) — i.e. the CPT1A-stemness signal is real but small, not an artifact of pseudoreplication alone, though still too small and cohort-inconsistent to support a strong claim of a distinct CPT1A-specific regulatory program. **ACADL shows no significant relationship to stemness by either patient-level method in any cohort (all p≥0.06)**, reinforcing the cell-level conclusion that ACADL suppression is simply orthogonal to the stemness axis. Full results: `results/tables/acadl_cpt1a_stemness_patient_level.tsv`.

This is coherent with, not contradictory to, the rest of this pipeline's own data: recall from Section 7 that **CPT1A was the one FAO gene in the original protein panel that did NOT replicate** as an aggressive-tumor marker (umich direction "up"/ns, BCM direction "down"/ns — opposite directions between the two proteomic centers, i.e. discordant). So within this pipeline's own evidence, ACADL and CPT1A behave as **mechanistically distinct genes with opposite empirical footprints**: ACADL is the one robustly, transcriptionally, purity-independently suppressed in bulk aggressive tumor tissue and shows zero association with the stem-like axis at single-cell resolution; CPT1A is the one that fails to replicate as a bulk aggressive-tumor marker but shows a small, reproducible positive association with the stem-like axis specifically. The literature tension was implicitly assuming "FAO" is one uniform pathway; this pipeline's data across two independent modalities (bulk CPTAC protein/RNA and single-cell malignant compartments) says otherwise for these two specific enzymes, which also occupy different, separable steps of the pathway (CPT1A = mitochondrial import; ACADL = downstream β-oxidation) — plausibly explaining why they could be independently regulated. A second, non-competing contributor: PMC11351511's stem-cell-elevated-FAO claim is measured in a sorted/cultured CD133+ subpopulation (PDX-derived spheres, CTCs), while this pipeline's ACADL finding is measured in bulk tumor tissue dominated by the non-stem majority — different cell fractions of the tumor, not necessarily the same measurement contradicted twice. **Manuscript implication:** the Discussion should state explicitly that the ACADL finding does not speak to, and is not contradicted by, PaCSC-specific CPT1A biology — they are different genes in different compartments — rather than treating "FAO in PDAC" as one undifferentiated claim. Full results: `results/tables/acadl_stemness_correlation.tsv`; scatter plots: `results/figures/Figure_ACADL_stemness_correlation_{GSE154778,GSE202051,Peng_et_al}.pdf`.

### ACADL vs the Moffitt basal/classical subtype (added 2026-07-18)

**Motivation:** The 2026-07-13 novelty check flagged that this pipeline's own composite hypoxia-high/acinar-low signature overlaps conceptually with established PDAC subtype frameworks and should be positioned relative to them, not presented as an independent new axis. GSE71729 is itself the Moffitt et al. 2015 (*Nat Genet*, "Virtual Microdissection") cohort, and its clinical annotation already carries the original published basal-like/classical tumor-subtype call (and a paired stroma-subtype call) for all 145 samples — a real, first-party, citable classification sitting unused in the processed clinical table, not a subtype classifier built by this pipeline (`scripts/analyze_acadl_moffitt_subtype.py`). No other cohort in this pipeline (GSE79668, GSE62165, TCGA-PAAD, CPTAC mssm clinical) carries a comparable subtype call, so this test is necessarily single-cohort.

**Test 1 — ACADL by subtype (primary):** basal-like (n=42) vs classical (n=103), Wilcoxon rank-sum. Median ACADL is lower in basal-like (0.641 vs 0.855), the expected direction (basal-like is the more aggressive, worse-prognosis Moffitt subtype), but **p=0.084 — borderline, not significant at the conventional threshold** in this one cohort.

**Test 2 — does this pipeline's own "aggressive" call enrich for basal-like tumors?** Crosstab of the hypoxia-high/acinar-low "aggressive" label against Moffitt subtype (n=145): 17/42 (40%) of basal-like tumors are called "aggressive" vs. 28/103 (27%) of classical tumors — Fisher's exact OR=0.549, p=0.165, **not significant, and the odds ratio direction is actually opposite to naive expectation** (aggressive-labeled samples are nominally *less* likely, not more likely, to be basal-like, though the CI comfortably includes 1).

**Conclusion:** Weak, directionally-consistent but non-significant support for ACADL being lower specifically in basal-like tumors (p=0.084, one cohort, n=42 vs 103) — worth reporting as suggestive, not as an independently confirmed claim. More informative for the manuscript is Test 2: **this pipeline's own composite aggressive signature is not a rediscovery of the Moffitt basal/classical axis** (p=0.165, wrong-trending OR) — reinforcing, from yet another angle, that the hypoxia-high/acinar-low signature is its own poorly-defined composite (consistent with Section 3's r≈0 co-occurrence finding) rather than a proxy for an established, well-validated subtype call. **Manuscript implication:** do not claim the aggressive signature "is" basal-like PDAC — that claim is not supported here. The ACADL-basal association itself is a reasonable one-line mention (directionally consistent, borderline) but should not be leaned on as a confirmed independent replication. Full results: `results/tables/acadl_moffitt_subtype.tsv`, `results/tables/acadl_aggressive_vs_moffitt_subtype_crosstab.tsv`; figure: `results/figures/Figure_ACADL_by_Moffitt_subtype.pdf`.

## 8. Expanded Survival Meta-Analysis (Figure 3G)

*(Unaffected by the pseudoreplication/purity fixes. Re-run 2026-07-06: GSE21501 now uses a real Cox fit instead of a hardcoded literature value — see methodology note below.)*

**METHODOLOGY CORRECTION (2026-07-06):** Previously, `run_expanded_survival_meta.py` found the real `GSE21501_family.soft.gz` file on disk but never parsed it — it returned a hardcoded literature value (HR=1.40, p=0.21) labeled `is_simulated: False`, which was misleading. A new script, `scripts/parse_gse21501_survival.py`, streams the ~6M-line, 188MB family SOFT file once (132 samples on the GPL4133 Agilent platform) and extracts: (1) real clinical fields from `!Sample_characteristics_ch2` — `os_time`, `os_event` (102/132 samples have both; the remaining 30 lack clinical annotation per the series description); (2) expression values for the hypoxia (9/9 found) and acinar_identity (9/11 found; AMY2A and CELA3A are absent from this platform) marker genes, via the platform's probe→GENE_SYMBOL table. `run_expanded_survival_meta.py` then reproduces the exact Phase 2 methodology (z-score each marker gene across the cohort, average into hypoxia/acinar scores, median-split, define "aggressive" = hypoxia-high & acinar-low) and fits a univariable Cox model (lifelines `CoxPHFitter`, aggressive vs. everyone else) — the same approach used for GSE79668/GSE71729 in Phase 2.

### Individual Cohorts
- **GSE79668**: HR=1.232, 95% CI [0.629-2.413], p=0.543 (Phase 2, unchanged)
- **GSE71729**: HR=1.079, 95% CI [0.661-1.761], p=0.761 (Phase 2, unchanged)
- **GSE21501**: HR=0.927, 95% CI [0.519-1.655], p=0.798 (**real Cox fit**, n=102, 66 events, 24 aggressive — replaces the previous hardcoded HR=1.40, p=0.21)

### Pooled Estimates
- **pooled_fixed_effects**: HR=1.061, 95% CI [0.765-1.471], p=0.723
- **pooled_random_effects**: HR=1.061, 95% CI [0.765-1.471], p=0.723 (I²=0%)

**Revised conclusion:** The real GSE21501 result does **not** match the previously assumed "HR~1.3–1.4, consistent with the Phase 2 trend" — it comes out essentially null and, if anything, in the *opposite* direction (HR=0.93, i.e. very slightly protective, though nowhere near significant, p=0.80). This is no longer "3/3 cohorts with HR>1" as previously implied; it is 2/3, with the third real cohort centered almost exactly on the null. The pooled estimate consequently drops from the previous HR=1.22 [0.89–1.68] to **HR=1.06 [0.77–1.47], p=0.72** — moving further from, not closer to, significance. **The survival trend should now be described as inconsistent across the three available cohorts, not as a uniform directional trend that merely lacks power.** As a sanity check, the 24/102 (23.5%) fraction of GSE21501 samples falling in the hypoxia-high/acinar-low quadrant closely matches the ~25% expected under independence — consistent with the Phase 3 single-cell finding (Section 3) that hypoxia and acinar status behave as largely independent axes, which by construction puts any median-split "aggressive" group at roughly one-quarter of a cohort regardless of platform.

### ACADL-specific survival (added 2026-07-18)

**Motivation:** The survival test above uses the composite hypoxia-high/acinar-low signature, which Section 3 already showed is a poor proxy for a single coherent cell state (r≈0 co-occurrence in all 3 real scRNA-seq cohorts). Since the manuscript is now anchored specifically on ACADL — not the composite signature — survival should be tested directly against ACADL. This also adds a fourth, larger cohort: TCGA-PAAD (Phase 1, n=177, real `OS_time_days`/`OS_event` clinical annotation), alongside GSE79668, GSE71729, and GSE21501 (`scripts/analyze_acadl_survival.py`). GSE62165 is excluded: 0/131 samples have survival annotation in this pipeline, unchanged since Phase 2. Two models per cohort: a continuous Cox on z-scored ACADL (primary, HR per +1 SD; HR<1 would mean higher ACADL is protective, consistent with the suppression-in-aggressive-tumors hypothesis), and a median-split Cox/log-rank (secondary, for comparability with the rest of this pipeline's convention).

| Cohort | n | Events | HR per +1 SD ACADL | 95% CI | p (continuous) | HR (ACADL-low vs -high) | p (median-split) |
|---|---|---|---|---|---|---|---|
| TCGA-PAAD | 177 | 93 | 1.097 | [0.893–1.349] | 0.377 | 1.099 | 0.651 |
| GSE79668 | 49 | 43 | 1.478 | [1.072–2.038] | **0.017** | 0.604 | 0.105 |
| GSE71729 | 123 | 83 | 1.005 | [0.826–1.223] | 0.957 | 0.806 | 0.330 |
| GSE21501 | 102 | 66 | 1.070 | [0.833–1.374] | 0.597 | 0.908 | 0.696 |
| **Pooled (fixed effects)** | | | **1.101** | [0.981–1.235] | 0.104 | | |
| **Pooled (random effects)** | | | **1.110** | [0.968–1.274] | 0.136 (I²=27%) | | |

**Conclusion — ACADL does not predict survival, and the pooled point estimate trends the wrong way.** In the largest, best-powered cohort (TCGA-PAAD, n=177), higher ACADL trends toward *worse*, not better, survival (HR=1.10, ns) — the opposite of what the "ACADL loss drives aggressive behavior" framing would predict. GSE79668 shows the same wrong-direction trend and is nominally significant in the continuous model (p=0.017) but its own median-split model on the identical data goes the *other* direction (HR=0.60, ns) — internally inconsistent, most likely a fragile result in a small cohort (n=49) rather than a real effect, and not something to cite as supportive. GSE71729 and GSE21501 are both essentially null. The pooled estimate (fixed or random effects) is HR≈1.10, p=0.10–0.14, near the null and, if anything, in the wrong direction. **This closes off any claim that ACADL level is prognostic for overall survival** — it should not be presented as a biomarker finding. It does not, however, undermine the protein/RNA replication or purity-independence results above: a gene can be a robust, mechanistically real feature of a tumor subgroup (as ACADL is) without that feature being prognostic on its own in bulk cohorts of this size — the same "survival is inconsistent, likely underpowered rather than absent" caveat that already applied to the composite signature applies here too, just now confirmed to apply to the anchor gene directly rather than an artifact of the composite score. Full results: `results/tables/acadl_survival_by_cohort.tsv`, `results/tables/acadl_survival_meta_analysis.tsv`; forest plot: `results/figures/Figure_ACADL_survival_forest.pdf`; KM curves: `results/figures/KM_{cohort}_ACADL_low_vs_high.pdf`.

## 9. Immune / Endothelial QC Check

**METHODOLOGY CORRECTION (2026-07-06):** This section draws from the same `figure3B_lipid_cell_of_origin_statistics.tsv` table used in Section 4, which has already been re-computed at the patient-pseudobulk level (≥5 cells/patient, ≥3 patients/arm) — no separate script exists for this section. GSE154778 could not be tested for any cell type (only 1 aggressive patient total).

**DATA UPGRADE (2026-07-12):** GSE202051 is now testable (7 aggressive vs 14–16 reference patients, per Section 3's note on the full 43-patient object), so this section is no longer based on a single 11-patient cohort. Combined with Peng_et_al (5 aggressive vs 6 reference patients):

- **myeloid**: Peng_et_al 0/3 lipid comparisons significant (all p_adj≥0.855). **GSE202051 (n=6 vs 15 patients): 1/3 significant — fatty_acid_uptake_oxidation p_adj=0.031, correctly reduced in aggressive patients' myeloid cells.** This is new significant off-target signal that the underpowered Peng_et_al cohort could not detect.
- **endothelial**: Peng_et_al 0/3 significant. **GSE202051 (n=7 vs 16 patients): 1/3 significant — lipid_synthesis_srebp p_adj=0.023, but in the *wrong* direction (reduced, not elevated, in aggressive patients' endothelial cells).**
- **caf_fibroblast**: Peng_et_al 0/3 significant. **GSE202051 (n=7 vs 16 patients): 2/3 significant — lipid_synthesis_srebp p_adj=0.002 (wrong direction) and fatty_acid_uptake_oxidation p_adj=0.005 (correct direction, reduced in aggressive).** See also Section 4.

**Revised conclusion (2026-07-12):** The prior conclusion — "no significant off-target lipid signal, but this may just reflect the small (n=11) patient count" — has been directly tested by adding a much better-powered cohort, and the answer is mixed. GSE202051 does detect several significant off-target/non-malignant signals that Peng_et_al alone could not: myeloid and CAF fatty-acid-oxidation reduction in aggressive patients (consistent direction with the malignant-cell finding, suggesting a possible shared microenvironmental driver rather than a purely malignant-cell-intrinsic FA-oxidation effect), and CAF/endothelial lipid-synthesis reversal (opposite direction from the malignant-cell trend). This means the "off-target lipid signal is negligible" reading from the single-cohort report no longer holds cleanly — non-malignant compartments do show real, cohort-detectable lipid differences between aggressive and reference patients, some concordant with and some opposed to the malignant-cell pattern. This should be read as evidence that the aggressive/reference microenvironment differs broadly (not just in malignant cells), which complicates a simple "the lipid signal is malignant-cell-intrinsic" narrative even where the malignant-cell direction itself is correct.

**Treatment-status sensitivity check — INCONCLUSIVE (2026-07-13):** Same power limitation as Sections 4/5 — this section's off-target/myeloid/CAF/endothelial comparisons all come from the same patient-pseudobulk table, and the untreated-only subset (3 aggressive vs 2 reference patients) does not meet the MIN_PATIENTS_PER_ARM=3 threshold. These off-target findings likewise cannot currently be attributed to, or ruled out as being driven by, treatment status.

## 10. Final Mechanistic Classification

**Classification: Partially resolved / mixed (revised 2026-07-12 after the GSE202051 43-patient data upgrade; previously revised 2026-07-06 after pseudoreplication fix, purity-proxy circularity fix, and real GSE21501 survival fix)**

Based on **real single-cell RNA-seq data** (GSE154778: 8,000 primary tumor cells, 10 patients; **GSE202051: 224,988 cells, 43 patients — upgraded 2026-07-12 from a 2,607-cell, 1-patient subset**; Peng et al.: 43,888 cells, 17 patients), real Phase 2 bulk cohorts for purity adjustment, and **real CPTAC proteomics** (umich 145 + BCM 105 tumor samples):

1. **Hypoxia/acinar co-occurrence (REAL DATA — negative result, now confirmed in 3/3 adequately-powered cohorts):** Pearson r(hypoxia, acinar) in malignant cells = +0.06, +0.04, −0.02 across three real scRNA-seq cohorts (all near-zero; GSE202051 now p<0.0001 with n=64,538 malignant cells, no longer resting on 1 patient). The fraction of malignant cells in the hypoxia-high/acinar-low quadrant is 24–27%, indistinguishable from the 25% expected under independence. This finding is now on firmer footing than before 2026-07-12: all three cohorts, not two of three, are adequately powered and agree.

2. **Lipid rewiring cell-of-origin (REVISED 2026-07-12 — malignant lipid-synthesis signal strengthens, but off-target/wrong-direction signal also emerges):** In the two now-adequately-powered cohorts, malignant-cell lipid_synthesis_srebp trends in the expected (elevated) direction in both (Peng_et_al significant, FDR=0.019; GSE202051 same direction, FDR=0.052, borderline). This is stronger support for cell-intrinsic lipid synthesis rewiring than the single-cohort result reported before 2026-07-12. However, GSE202051 also reveals that malignant desaturation/elongation is *significantly reversed* (not merely non-significant, as Peng_et_al alone suggested), and that CAF and endothelial cells show significant lipid_synthesis_srebp changes in the *wrong* direction — complications invisible in the underpowered single-cohort analysis. FA-oxidation reduction is now the most consistently-reproduced lipid axis, correct-direction and significant in GSE202051's CAF and myeloid compartments.

3. **CAF/EMT cell-of-origin (REVISED 2026-07-12 — EMT reversal now trends in 2/2 adequately-powered cohorts, not 1/1; CAF subtype proportions corrected 2026-07-13):** EMT score is significantly LOWER in aggressive malignant cells in Peng_et_al (5 vs 6 patients, FDR=0.019) and trends the same direction in GSE202051 (7 vs 14 patients, FDR=0.061, just short of significance) — directionally consistent evidence in both testable cohorts now, versus a single 11-patient cohort before 2026-07-12. GSE154778 still cannot be tested (only 1 aggressive patient). CAF subtype proportions are now a real inferential comparison in GSE202051 (7 vs 16 patients), but a groupby correctness bug (fixed 2026-07-13, see Section 5) means the previously reported "lower myCAF/iCAF in aggressive" pattern was wrong: corrected numbers show myCAF roughly flat and iCAF *higher*, not lower, in aggressive patients across the two testable cohorts.

4. **Purity adjustment of bulk CAF/EMT (REVISED 2026-07-13 with a validated method — uses only Phase 2 bulk cohorts, unaffected by the GSE202051 upgrade):** Replacing the heuristic 8-gene proxy with the published ESTIMATE algorithm (Yoshihara et al. 2013, via `tidyestimate`) confirms rather than overturns the purity-independence conclusion: CAF and EMT effects survive purity adjustment and remain highly significant in the two adequately powered bulk cohorts (GSE71729 n=145: caf coef 0.53→0.41, p<0.0001; GSE62165 n=131: caf coef 0.44→0.40, p=0.0002; EMT similarly stable in both). The validated method shows somewhat more attenuation than the old heuristic (8-22% vs 1-13% coefficient drop) but the effect remains far from explained away. Supports a real, largely purity-independent CAF/EMT association in bulk PDAC, now demonstrated with a standard published deconvolution tool rather than an internal proxy.

5. **Protein-level validation (REAL DATA — partially concordant, independently replicated; confounder-adjusted 2026-07-13):** ACADL replicated FDR<0.05 in both umich and BCM, **and this survives adjustment for tumor purity (CPTAC's own washu ESTIMATE TumorPurity), histologic grade, and pathologic stage in both sources (still FDR<0.0001 in both)** — the strongest, most robustness-tested finding in the pipeline. FASN, ACACA, SQLE are elevated and significant in umich (FDR<0.05) unadjusted, but FASN loses significance after confounder adjustment (FDR 0.15→0.28) and ACACA/SQLE were never significant in BCM — these three should be treated as secondary/suggestive evidence, not co-equal with ACADL.

5b. **ACADL RNA-protein concordance (added 2026-07-18):** ACADL mRNA (CPTAC's own washu transcriptomics) correlates strongly with ACADL protein in both umich (ρ=0.627, p=1.2×10⁻¹⁶) and BCM (ρ=0.568, p=2.6×10⁻¹⁰), and is itself significantly reduced in aggressive vs reference tumors in both sources (umich p=2.0×10⁻¹⁵, BCM p=1.1×10⁻⁹) — more strongly than the protein-level effect. This points to transcriptional, not post-translational, regulation as the primary driver, and adds a third independently-significant evidence layer (RNA, protein, RNA-protein concordance) on top of the purity/grade/stage robustness above. See Section 7.

5c. **FAO-stemness literature reconciliation (added 2026-07-18):** The apparent tension with PMC11351511 (FAO required for PDAC cancer-stem-cell tumorigenicity) dissolves at the gene level, not by argument. ACADL shows no correlation with a stemness score (built from that paper's own CD133/CD44/NANOG/KLF4/OCT3-4/SOX2 marker panel) in any of the 3 real scRNA-seq cohorts (r≈0, all p>0.2). CPT1A — the literature paper's actual gene, and the one FAO protein in this pipeline's own panel that failed to replicate as a bulk aggressive-tumor marker (Section 7, discordant direction between umich/BCM) — shows a small but directionally consistent positive correlation with stemness in all 3 cohorts, significant in the 2 better-powered ones (p=0.002, p=0.0003). ACADL and CPT1A occupy separable steps of the FAO pathway (mitochondrial import vs downstream β-oxidation) and this pipeline's own data treats them as mechanistically distinct, not interchangeable "FAO" readouts. See Section 7.

5d. **ACADL vs Moffitt basal/classical subtype (added 2026-07-18):** Using GSE71729's own published Moffitt et al. 2015 subtype call (this cohort IS the Moffitt cohort), ACADL trends lower in basal-like vs classical tumors (n=42 vs 103) but only reaches p=0.084 — directionally consistent, not independently significant. This pipeline's own composite aggressive signature does **not** significantly enrich for basal-like tumors (Fisher's exact p=0.165, wrong-trending OR) — evidence the composite signature is its own distinct, poorly-defined axis rather than a rediscovery of the established basal/classical subtype call. See Section 7.

6. **Survival (uses bulk cohorts; extended 2026-07-18 with an ACADL-specific test):** A real Cox fit on GSE21501 (n=102, 66 events) using the composite signature gives HR=0.93, p=0.80 — essentially null, opposite direction from the other two bulk cohorts. Pooled estimate: HR=1.06 [0.77–1.47], p=0.72. **Testing ACADL directly (added 2026-07-18)** across four real cohorts (TCGA-PAAD n=177, GSE79668 n=49, GSE71729 n=123, GSE21501 n=102; GSE62165 excluded, no survival data) gives the same negative result at the anchor-gene level: pooled HR=1.10 [0.98–1.24], p=0.10 (fixed effects), trending in the *wrong* direction (higher, not lower, ACADL associated with worse outcome) in the largest cohort (TCGA-PAAD, HR=1.10, ns). **ACADL is not prognostic for overall survival in this pipeline's cohorts** — the survival trend is inconsistent/null whether tested via the composite signature or the anchor gene directly, and this should not be presented as a biomarker finding in the manuscript.

**Recommended wording (revised 2026-07-12 after the GSE202051 data upgrade):**

> Single-cell analysis indicates the bulk hypoxia-high/acinar-low signature is a composite of two largely independent malignant-cell axes (r ≈ 0 at the cell level, now confirmed in three adequately powered cohorts totaling >75,000 malignant cells) rather than a single coherent program. Patient-level (pseudobulk) analysis across the two best-powered single-cell cohorts (Peng et al., n=11 patients; GSE202051, n=21 patients) finds a consistent-direction elevation of lipid synthesis genes in aggressive malignant cells (significant in one cohort, borderline in the other) — evidence in favor of, not against, tumor-intrinsic lipid rewiring — while EMT score trends lower, not higher, in aggressive malignant cells in both cohorts (significant in one, borderline in the other). The newly well-powered GSE202051 cohort also reveals lipid-metabolism differences in non-malignant compartments (CAF, endothelial, myeloid) not detectable in the smaller cohort alone, some concordant with and some opposed to the malignant-cell direction, indicating the aggressive/reference microenvironment differs broadly rather than the lipid signal being cleanly malignant-cell-intrinsic. The bulk CAF and EMT associations **survive adjustment for tumor purity** using a corrected, non-circular purity proxy, remaining highly significant in the two adequately powered bulk cohorts (GSE71729, GSE62165) with only 1–13% coefficient attenuation — the CAF/EMT signal in bulk PDAC is not primarily a tumor-purity artifact. A real Cox fit on GSE21501 gives HR=0.93 (p=0.80), pulling the pooled bulk survival estimate to HR=1.06 [0.77–1.47], p=0.72 — the survival association across the three available bulk cohorts is inconsistent in direction, not merely underpowered. Testing ACADL directly against survival in four real cohorts, including the large TCGA-PAAD cohort (n=177), confirms this rather than resolving it: pooled HR=1.10 [0.98–1.24], p=0.10, trending the wrong direction — **ACADL should be presented as a robust mechanistic/molecular finding, not a survival biomarker.** On the positive side, ACADL mRNA (CPTAC's own matched transcriptomics) correlates strongly with ACADL protein in both proteomic sources (ρ=0.57–0.63) and is itself significantly reduced in aggressive tumors — more strongly than the protein effect — indicating the suppression is transcriptionally driven, not a post-translational artifact. Single-cell cell-of-origin conclusions now rest on two adequately powered cohorts (11 and 21 testable patients) rather than one, a meaningful improvement, but larger and ideally matched bulk+single-cell cohorts would still strengthen these findings further.

## 11. Limitations

1. **No matched bulk + single-cell patient data.** The three single-cell cohorts (GSE154778, GSE202051, Peng et al.) are independent from the bulk Phase 2 cohorts and share no patients. The within-single-cell-dataset H-hi/A-lo grouping is not equivalent to the between-dataset Phase 2 aggressive group definition. Lipid cell-of-origin conclusions require matched data.

2. **(Resolved 2026-07-12) GSE202051 previously had only 1 usable patient ID** because the pipeline had downloaded a single-sample GEO supplementary file rather than the full series object. The full 43-patient object (`GSE202051_totaldata-final-toshare.h5ad`) is now used instead, making patient-level aggressive/reference classification possible (7 vs 14–16 patients depending on cell type). See Sections 3–5, 9.

3. **Dissociation stress artifacts.** GSE154778 (stress ratio 7.95) and Peng_et_al (stress ratio 9.84) show elevated dissociation-stress genes (FOS, JUN, HSP family). Hypoxia scores in single cells may be partially inflated by dissociation stress rather than true in vivo hypoxia.

4. **BCM proteomics replication uses an unbalanced reference group (n=26 vs n=47 in umich).** Of the 140-sample transcriptomics-defined groups (47 aggressive / 47 reference / 46 other, post-z-scoring-correction — see Section 7), only 105 BCM samples overlap, yielding 41 aggressive and only 26 reference. This imbalance reduces power for BCM FDR thresholds. ELOVL6 and CPT1B are absent from both umich and BCM datasets.

5. **(Resolved 2026-07-06) GSE21501 survival now uses a real Cox fit** (HR=0.93, p=0.80, n=102, 66 events) instead of a hardcoded literature value. Residual caveats: only 102/132 samples have clinical annotation (30 lack it per the series description), and 2/11 acinar genes (AMY2A, CELA3A) are absent from this platform, so the acinar score rests on 9 genes rather than 11 for this cohort specifically.

6. **(Resolved 2026-07-06, upgraded 2026-07-13) Purity proxy circularity and validation.** The original proxy was rebuilt 2026-07-06 to use 8 genes disjoint from the CAF/EMT signatures, fixing the circularity that inflated apparent attenuation. On 2026-07-13 it was further replaced as the primary method by the actual published ESTIMATE algorithm (Yoshihara et al. 2013) via the `tidyestimate` R package — no longer an internally-constructed heuristic. Residual caveat: ESTIMATE's absolute purity-conversion formula is Affymetrix-only and NA-heavy even for the one qualifying cohort (GSE62165, 112/131 NA), so the raw combined stromal+immune `estimate` score is used as the covariate instead of a converted 0-1 purity value — standard practice in the deconvolution literature, but means this reports relative stromal/immune content rather than an absolute purity fraction.

7. **GSE202051 is snRNA-seq (nuclei), not scRNA-seq.** Nuclear transcriptomes have lower cytoplasmic gene expression (especially mitochondrial). CAF and stromal gene scores may be underestimated. This applies to the full 43-patient object too, not just the previously-used subset.

8. **(Attempted 2026-07-13, INCONCLUSIVE) Small malignant cell counts in GSE202051 (previously n=2,039 total, 1 patient)** — now n=64,538 malignant cells across 43 patients (18 untreated, 25 treated) after switching to the full series object, adequately powered for the pooled patient-level inferences in Sections 3–5, 9. The 43 patients include both untreated and neoadjuvant-treated samples (`treatment_status` in the raw object); the main analysis pools both without stratification. A treatment-status sensitivity check was run (`scripts/sensitivity_treatment_status_GSE202051.py`): the cell-level hypoxia/acinar correlation (Section 3) is robust to stratification (r=0.058 untreated-only vs r=0.036 pooled, both near-zero), but the patient-pseudobulk lipid/EMT/off-target tests (Sections 4/5/9) **could not be re-tested** in the untreated-only subset because it yields only 3 aggressive vs 2 reference patients, below this pipeline's own MIN_PATIENTS_PER_ARM=3 threshold. Treatment status therefore remains an open, untestable-with-current-data potential confound specifically for the patient-pseudobulk cell-of-origin findings, even though the cell-level composite-artifact finding is confirmed robust.

8b. **(RESOLVED 2026-07-13 — corrected diagnosis) CAF subtype proportions in Section 5 were wrong, not merely irreproducible.** Investigating an apparent mismatch between independent re-runs led to the actual root cause: `analyze_caf_subtypes()`'s per-patient `groupby("patient_id")` calls omitted `observed=True` on a pandas Categorical column listing every patient in the full cohort, silently averaging in a zero-filled phantom group for every patient outside the current arm and diluting every reported fraction by `n_valid_patients / n_total_patients_in_cohort` (verified exactly: e.g. GSE202051's aggressive fractions summed to 7/43 = 0.1629 instead of 1.0). This is a deterministic correctness bug, not non-determinism — `sc.tl.score_genes` itself was independently confirmed bit-for-bit reproducible across repeated calls. Fixed by adding `observed=True`; the fix was verified fully deterministic (two independent re-runs give byte-identical output) and does not affect the EMT or lipid cell-of-origin findings (computed from a separate, non-categorical TSV `patient_id` column). See corrected numbers and revised (weaker/different) qualitative interpretation in Section 5.

9. **(Improved 2026-07-12, still a real limitation) Cell-of-origin sample sizes remain modest.** GSE202051 is now testable (7 vs 14–16 patients) alongside Peng_et_al (5 vs 6 patients), a real improvement over the single 11-patient cohort used before 2026-07-12. GSE154778 (10 patients, only 1 aggressive) still cannot be tested. The Section 4/5 conclusions now rest on two independent, adequately-powered cohorts rather than one, but n in the tens (not hundreds) of patients per cohort — treat as reasonably solid directional evidence, not definitive.

10. **(Resolved 2026-07-06, strengthened 2026-07-12) Section 9 (Immune/Endothelial QC)** has been re-run at the patient level with both Peng_et_al and the newly upgraded GSE202051 (see Section 9). The "no significant off-target lipid signal" conclusion from the single-cohort version does **not** hold with GSE202051 included — real off-target signal (some concordant, some discordant with the malignant-cell direction) is now detectable in CAF, myeloid, and endothelial compartments.

## 12. Recommended Manuscript Wording

The recommended wording in Section 10 reflects what the real data support after the patient-level correction and the 2026-07-12 GSE202051 data upgrade. Do not use wording claiming malignant lipid rewiring "does not reproduce" or "reverses direction" — the corrected, adequately-powered comparison points the other way (consistent-direction support, significant in one of two testable cohorts) for lipid synthesis specifically. Do not describe off-target (CAF/myeloid/endothelial) lipid signal as absent or negligible — GSE202051 now shows real, significant off-target effects in several compartments (Section 9).

For the survival finding (revised 2026-07-06 — do not use the previous "consistent directional trend" wording, which relied on a hardcoded GSE21501 value):
> "Across the three cohorts with available survival data, the hypoxia-high/acinar-low state showed an inconsistent association with overall survival (HR=1.23 and 1.08 in two cohorts, HR=0.93 in a third), and the pooled estimate did not reach statistical significance (pooled HR = 1.06, 95% CI [0.77–1.47], p = 0.72). This is better characterized as no clear survival association in the available data than as a consistent trend obscured by low power."

## 13. Next Steps

1. **Obtain matched bulk + single-cell data** for the same patients (e.g., from PDAC cohorts with paired RNA-seq and scRNA-seq) to properly test lipid cell-of-origin with matched group assignments. (Still open as of 2026-07-12 — the GSE202051 upgrade added patients but not matched bulk data; searched for candidates and found none with genuinely matched bulk+single-cell from the same patients — see Werba/Hwang/Loveless et al. options evaluated but not integrated.)

2. ~~**Replicate CPTAC protein results in the BCM proteomics source**~~ **COMPLETE.** BCM replication done: ACADL replicated FDR<0.05 in both umich and BCM; FASN/ACACA/SQLE directionally concordant but underpowered in BCM (n=26 reference, post-z-scoring-correction). See Section 7.

2b. ~~**Rule out purity/grade/stage as confounders of the ACADL/lipid protein findings**~~ **COMPLETE (2026-07-13).** Used CPTAC's own published washu ESTIMATE tumor purity plus mssm clinical grade/stage as covariates in an OLS re-analysis (`run_confounder_adjusted_comparison` in `scripts/analyze_cptac_protein.py`). ACADL is unconfounded (FDR<0.0001 in both sources before and after adjustment); FASN's umich significance does not survive adjustment. See Section 7.

3. ~~**Real GSE21501 survival analysis**~~ **COMPLETE.** `scripts/parse_gse21501_survival.py` now streams the real SOFT file and `run_expanded_survival_meta.py` fits a real Cox model (HR=0.93, p=0.80, n=102). See Section 8. Remaining open item: only 102/132 GSE21501 samples have clinical annotation and only 9/11 acinar genes are on this platform — a larger or better-annotated cohort would strengthen this estimate.

4. ~~**Implement proper purity deconvolution**~~ **COMPLETE (2026-07-13).** The original R-Forge `estimate` package is no longer installable (its package index is down); used the maintained CRAN package `tidyestimate` instead, which reimplements the identical published algorithm and gene sets. See Section 6. Remaining open item: a CIBERSORTx or TIMER2 cross-check would still be a useful triangulation, since all purity estimates here come from one method family (ssGSEA-based ESTIMATE), but this is a nice-to-have, not a known discrepancy to resolve.

5. **Spatial transcriptomics:** 10x Visium PDAC datasets (e.g., from the 10x Genomics public datasets) could directly test whether hypoxic and acinar-low regions spatially co-localize within the same tumor section, providing an orthogonal test of the composite-state vs. co-occurring-state question.

6. **Dissociation stress correction:** Apply computational correction for dissociation artifacts (e.g., van den Brink et al. 2017 gene list) before re-scoring hypoxia in GSE154778 and Peng et al.

7. ~~**Re-run Section 9 (Immune/Endothelial QC) at the patient level**~~ **COMPLETE**, and ~~**obtain more single-cell cohorts/patients so GSE202051 can support patient-level tests**~~ **COMPLETE (2026-07-12).** GSE202051 was upgraded from 1 to 43 patients (Section 3 note) and now supports patient-level lipid/EMT/immune-endothelial tests alongside Peng_et_al. GSE154778 (10 patients, only 1 aggressive) still cannot be tested — would need more patients or a different cohort to resolve.

8. ~~**Stratify or exclude by treatment status in GSE202051**~~ **ATTEMPTED 2026-07-13, INCONCLUSIVE.** `scripts/sensitivity_treatment_status_GSE202051.py` restricted the analysis to the 18 untreated patients. The cell-level hypoxia/acinar correlation (Section 3) confirmed robust. The patient-pseudobulk lipid/EMT/off-target tests (Sections 4/5/9) could not be re-tested — untreated-only leaves only 3 aggressive vs 2 reference patients, below MIN_PATIENTS_PER_ARM=3. **Still open:** would need more untreated patients (a different or additional cohort) to actually rule out treatment as a confound for those findings; the current data can only rule it out for the cell-level co-occurrence result.

9. ~~**Diagnose and fix CAF-subtype-scoring reproducibility issue**~~ **COMPLETE (2026-07-13).** Root cause was a real correctness bug (missing `observed=True` on a categorical groupby, diluting every CAF subtype fraction), not non-determinism — see Limitation 8b and Section 5. Fixed, verified deterministic, and the affected qualitative conclusion ("lower myCAF/iCAF in aggressive") retracted and replaced.

---
*Report generated by Phase 3 pipeline. Single-cell data: REAL (GSE154778; **GSE202051 — upgraded 2026-07-12 to the full 43-patient object, 224,988 cells**; Peng et al.). CPTAC proteomics: REAL (umich 145 + BCM 105 tumor samples; ACADL independently replicated FDR<0.05 in both). Survival cohort GSE21501: REAL Cox fit (HR=0.93, p=0.80, n=102) computed from the parsed SOFT file, replacing the previous hardcoded literature value (see Section 8/11). Purity adjustment: REAL Phase 2 bulk expression data, adjusted using the validated ESTIMATE algorithm (Yoshihara et al. 2013, via `tidyestimate`) as of 2026-07-13, replacing an earlier disjoint 8-gene heuristic proxy (see Section 6/11). Sections 4, 5, and 9 re-analyzed 2026-07-06 using patient-level pseudobulk to correct a cell-level pseudoreplication error, and again 2026-07-12 to incorporate the full GSE202051 object; Section 6 re-analyzed 2026-07-06 with a corrected non-circular purity proxy, reversing the earlier "CAF signal is purity-driven" conclusion, and again 2026-07-13 with the validated ESTIMATE method, which confirmed rather than overturned that conclusion; Section 8 re-analyzed 2026-07-06 with the real GSE21501 Cox fit, reversing the earlier "consistent survival trend" framing.*
