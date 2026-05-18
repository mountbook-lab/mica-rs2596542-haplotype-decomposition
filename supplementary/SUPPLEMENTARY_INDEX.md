# Supplementary materials index — pre-AJHG submission

Date: 2026-05-15.  Project: `mica-rs2596542-haplotype-decomposition`.

All numbers in the captions below are taken directly from the corresponding table CSV. Source scripts are in the parent repo.

---

## Supplementary Figures (7)

**Figure S1 — Comparison of NMF and LDA decomposition on the same rs2596542-T carrier haplotype matrix.**
Source: `supp_1_lda_comparison.py`. k = 8 LatentDirichletAllocation fitted on the same 2,111 × 7,116 haplotype matrix; Hungarian matching of LDA components to NMF components via H-row Pearson correlation. Panel (a) — matched correlation per NMF component; the dashed red line marks the matching threshold r = 0.7. Panel (b) — LDA-c5 analogue eQTL / C-flip enrichment alongside NMF-c5 (c5-like enrichment pattern is *partially retained* under LDA); cells annotated `n/a` indicate zero overlap (fold = 0) between the matched component's top-5 % SNV set and the gene's eQTL set in the universe (no MICA eQTL overlap in either method's c5-like component). **Key values**: matched corr c4 / c5 / c6 = 0.99 / 0.87 / 0.96; LDA-c5 preserves MICA depletion (fold = 0) and C-flip enrichment (fold = 1.39) but attenuates the localized HLA-B *cis*-eQTL signal (NMF fold 2.14 → LDA fold 0.66). LDA recovered the major component structure but smoothed the localized HLA-B regulatory signature, supporting NMF as the appropriate method for sparse, localized component-defining SNV signatures.

**Figure S2 — NMF rank selection by consensus cophenetic correlation.**
Source: `supp_2_cophenetic_k.py`. k = 2..12 with 20 random initializations per rank; consensus matrix C[i,j] = (fraction of seeds where haplotypes i and j share dominant component); cophenetic correlation = Pearson r between cophenetic distances on the average-linkage dendrogram and (1 − C). Panel (a) — cophenetic correlation and Brunet-2004 dispersion vs k; (b) — mean reconstruction error ‖X − WH‖_F. **Key values**: cophenetic max = 0.998 at k = 2; k = 8 cophenetic = 0.963 (within 0.05 of max).

**Figure S3 — NMF seed reproducibility at k = 8.**
Source: `supp_3_seed_stability.py`. 50 independent NMF runs at k = 8 (init = random). Hungarian-matched H-row Pearson correlations (anchored to seed 0); recovery of MICA-enriched (fold > 3) and HLA-B-enriched (fold > 1.5) component labels per seed. Panel (a) — boxplot of matched corr per component; (b) — MICA / HLA-B fold across 49 matched seeds. **Key values**: c4 / c5 / c6 mean matched corr = 0.98 / 0.96 / 0.99; MICA-enriched recovered in 46/49 (94 %) seeds; HLA-B-enriched recovered in 49/49 (100 %) seeds.

**Figure S4 — LIRI-JP c5 score weighting sensitivity.**
Source: `supp_4_liri_weighted.py`. Three c5 score definitions are compared: unweighted top-5 % SNV burden (n = 335 NMF SNVs available in LIRI), H-loading-weighted across all 7,116 NMF SNVs, and NNLS projection of LIRI dosage onto the 1000G NMF basis. Panel (a) — unweighted vs weighted scatter; (b) unweighted vs NNLS; (c) forest plot of c5 β on log₂(HLA-C TPM) under the full M3 model (c5 + clinical + immune_PC1 + log₂(HLA-A); the same model as Table S1 row B and Supplementary Figure S5). **Key values** (n = 117; complete-case from 122 LIRI-JP donors with matched genotype + RNA-seq, after dropping 5 donors with missing T_stage or non-standard viral_status): Pearson r between unweighted and weighted = 0.79; β (unweighted / weighted / NNLS) = −0.193 / −0.160 / −0.175; all P < 0.005 in M3.

**Figure S5 — HLA-A adjustment diagnostic for the LIRI-JP regression of HLA-C on c5.**
Source: `supp_5_liri_hlaa_conditioning.py`. Diagnostic for whether HLA-A acts as a confounder of the c5 → HLA-C ↓ association or as a precision-improving covariate capturing pan-class-I shared expression variance. Panel (a) — c5 score vs log₂(HLA-A); (b) — log₂(HLA-A) vs log₂(HLA-C); (c) — c5 score vs HLA-C residual after regressing out HLA-A; (d) — forest plot of c5 β across M0–M3. **Key values** (n = 117; complete-case from 122 LIRI-JP donors with matched genotype + RNA-seq): c5 ⊥ HLA-A (Pearson r = −0.04, p = 0.64); HLA-A explains 70 % of HLA-C variance; c5 vs HLA-C residual r = −0.37, p = 5 × 10⁻⁵ (stronger than the raw r = −0.22); SE shrinkage M2 → M3 = −30 %; β change = +19 % (suppression direction); VIF(c5) = 1.10, VIF(HLA-A) = 2.10.

**Figure S6 — LIRI-JP HLA-C regression sensitivity to analytic-sample choice and multiple-testing context.**
Source: `supp_6_liri_sensitivity.py`. Three sensitivities: S1 (whitespace-trim of `viral_status` strings, n: 117 → 120), S2 (S1 + re-assigning HBV+HCV co-infected donors to HCV, n → 122), and S3 (c5 × viral_status interaction term on the primary n = 117 set). Panel (a) — sensitivity forest plot; (b) — −log₁₀(p) for the 7-outcome family (Table S1 section C) under raw, BH-FDR, and Bonferroni adjustment. **Key values**: primary β = −0.191 (p = 5×10⁻⁴), S1 β = −0.185 (p = 7×10⁻⁴), S2 β = −0.187 (p = 5×10⁻⁴); c5 × HCV interaction term β = +0.18, p = 0.24 (NS); c5 × NBNC term β = −0.05, p = 0.78 (NS). HLA-C is the only outcome that survives both BH-FDR and Bonferroni (adjusted p = 0.0037 by both methods, α = 0.05); the six other tested gene outcomes are pure null after correction. Locus specificity confirmed.

**Figure S7 — Reference-panel replication and NMF extension of the COMT Val158Met flip-flop example.**
Source: `figure_S7_comt_flipflop.py` (consumes `~/analyses/comt_flipflop_nmf/tables/comt_pair_metrics_by_pop.tsv`, `comt_component_mixture_by_rs4680_status.tsv`, `comt_component_mixture_by_pop.tsv`). (A) Population-conditional signed LD and asymmetry component C for rs4680 (Val158Met) and rs2097603 across 1000 Genomes populations. Blue and brown indicate positive and negative signed values, respectively. (B) Row-normalized k = 4 NMF component mixture conditional on rs4680 allele status, with V indicating Val/REF carriers and M indicating Met/ALT carriers. (C) Population-level COMT NMF component mixture across all 5,008 haplotypes. Compared with the MHC rs2596542 locus, the COMT decomposition is simpler and more population-like, supporting broader applicability of carrier-set decomposition while illustrating locus-dependent differences in functional interpretability. COMT is included as an external reference-panel proof-of-concept and is not cited in the main Results.

## Supplementary Tables (10)

**Table S1 — LIRI-JP regression models.**
Built by `make_supplementary_package.py` from `liri_c5_score.csv`, `EGA_clinical_matched.csv`, and `tpm_gene_named.csv`. Multi-section table:

- *Section A*: HLA-C ~ c5 (unweighted) progressive adjustment M0 (c5 only) → M1 (+ clinical) → M2 (+ immune_PC1) → M3 (+ HLA-A).
- *Section B*: HLA-C M3 full model with three c5 score variants (unweighted / H-loading-weighted / NNLS-projected).
- *Section C*: M3 full model fit to each MHC class-I / class-III outcome (HLA-B, HLA-C, MICA, MICB, HCG27, HCP5, HLA-A).
- Companion `table_S1d_VIF.csv` reports VIF for every predictor in the full M3 design matrix (c5 + clinical + immune_PC1 + HLA-A).
  **Rows**: 14 regression rows + 8 VIF rows.

**Table S2 — Top-signature SNVs defining NMF components in rs2596542-T carrier haplotypes.**
Built by `make_supplementary_package.py` from `results/nmf_H_26_k8.parquet` and `results/region_per_partner_26.csv`. Top-signature SNVs were defined as variants with H-matrix loadings in the top 5 % of each component's H[c, ·] distribution (i.e., at or above the 95th percentile), yielding 355 SNVs per component (8 × 355 = 2,840 rows). The table lists the component-defining SNV sets used for LD-reversal enrichment, *cis*-eQTL annotation, MAF×distance matched-permutation analysis (Supplementary Table S3), and tag-SNP component mapping (Supplementary Table S8). The anchor SNV rs2596542 was removed from the NMF input matrix after carrier selection because all retained haplotypes carry the rs2596542-T allele at this site, leaving it with zero variance in the analytic matrix. Columns: `component` (c0..c7); `rank_within_component` (1 = highest H_loading, descending); `position_grch37` (chr6 1000G phase 3 coordinate); `H_loading` (rounded to 4 decimal places); `gene_region` (MICA / MICB / HLA-B / HLA-C / HCG27 / HCP5 / 3'_of_MICB / 5'_of_HCG27 / intergenic_class_I); `C_flip`, `class_flip` (boolean LD-reversal flags); `eQTL_<gene>` (boolean indicating whether the SNV is a GTEx Whole_Blood *cis*-eQTL for that gene). **Rows**: 2,840 (8 components × 355).

**Table S3 — Matched permutation for component-level enrichment.**
Source: `permutation_all_components.py`. Per (NMF component × test set) MAF×distance matched permutation, 10,000 perms. Test sets: r-flip / C-flip / class-flip / MICA / MICB / HLA-B / HLA-C eQTL. Columns: component, test_set, n_top, obs_overlap, matched_mean, matched_sd, fold_vs_matched, z_score, empirical_p (one-sided), empirical_p_reported (string). MAF bins [0, 0.01, 0.05, 0.10, 0.20, 0.30, 1.0]; distance bins [0, 25, 50, 100, 150, 200, 260] kb. *Legend note*: empirical *p* = 0 in the numeric column indicates that no permutation exceeded the observed statistic and is reported as **P < 1×10⁻⁴** (the 10,000-permutation floor) in the `empirical_p_reported` column.

**Table S4 — c5 axis eQTL direction across GTEx tissues.**
Source: `eqtl_axis_test.py`, `liver_replication.py`, `external_tissue_replication.py`. Per (analysis_block, component, tissue, target_gene): n_overlap, mean / median NES, n_pos / n_neg / sign_coherence (= |mean(sign(NES))|). Two analysis blocks (column `analysis_block`):

- *all_components_primary*: 8 NMF components × 5 target genes, fitted in **Whole_Blood** and **Liver**.
- *c5_cross_tissue_replication*: c5 only × 5 target genes, fitted in **Cells_EBV-transformed_lymphocytes**, **Lung**, and **Skin_Sun_Exposed_Lower_leg** to test reproducibility of the c5 direction in additional tissues.

  **Rows**: 87. The c5 rows in the cross-tissue block are *not* duplicates of the Whole_Blood / Liver primary block; they describe the same component evaluated in different tissues.

**Table S5 — LDA vs NMF c5-analogue enrichment.**
Source: `supp_1_lda_comparison.py`. Side-by-side fold enrichment for NMF c5 and LDA c5-analogue (Hungarian-matched) on 6 test sets (5 eQTL gene sets + C-flip).

**Table S6 — NMF rank cophenetic and reconstruction error.**
Source: `supp_2_cophenetic_k.py`. Per k ∈ {2, 3, ..., 12}: cophenetic_corr, dispersion, mean_recon_err, std_recon_err over 20 random-initialization seeds.

**Table S7 — NMF seed stability per component.**
Source: `supp_3_seed_stability.py`. Per component (c0..c7): mean / min / std of Hungarian-matched H-row Pearson correlation across 49 seeds vs seed 0 reference. Companion `table_S7b_nmf_seed_stability_per_seed.csv` lists per-seed matched correlations and MICA / HLA-B label recovery flags.

**Table S8 — Tag SNP → NMF component mapping.**
Source: `tag_snp_component_map.py`, with biological-role annotations added by `make_supplementary_package.py`. Four manuscript-relevant tag SNPs (`rs2596542` anchor; `rs2395029` HLA-B*57:01 surrogate in HCP5; `rs2244546` Lange-2013 HCV-HCC tag in HCP5; `rs11509487` MICB upstream indel — the paper-1 anchor) are intersected with the 8 NMF components. The `in_NMF_set` column is True for `rs2395029` and `rs2244546`, which are biallelic SNVs that pass the MAF ≥ 0.05 filter and are direct rows of the H matrix. It is False for `rs2596542` (anchor, excluded from H by construction; used only to define the carrier set) and for `rs11509487` (a 9-bp insertion filtered from the biallelic SNV matrix). For these two, the H values reported are not for the tag itself but for the *nearest in-window NMF SNV* (`nearest_NMF_SNV_position` column); the column `distance_bp_to_nearest_NMF_SNV` makes this explicit (rs2596542: 374 bp; rs11509487: 43 bp). Two CSV variants are provided:

- `table_S8_tag_snp_component_mapping.csv` — long format (1 row per (tag SNP × component); 32 rows). Columns: rsID, biological_role, grch37_pos, in_NMF_set, component, H_loading, rank_within_component_descending, rank_percentile_within_component, in_top_5pct_of_H_loadings, nearest_NMF_SNV_position, distance_bp_to_nearest_NMF_SNV, distance_bp_to_nearest_top5pct_SNV.
- `table_S8_tag_snp_component_mapping_wide.csv` — wide format (1 row per tag SNP; 4 rows). Columns: rsID, biological_role, grch37_pos, in_NMF_set, H_c0..H_c7, in_top5pct_c0..c7, best_component, best_H_loading, n_components_in_top5pct.

Manuscript-relevant findings (best component = argmax H_loading; for tag SNPs not in the NMF matrix, H values are at the nearest in-NMF SNV — see column `distance_bp_to_nearest_NMF_SNV`):

- **rs2596542** (anchor; excluded from H by construction; nearest NMF SNV 374 bp away): nearest SNV's argmax is c3 (H = 0.33, rank 13 % within c3); H_c5 = 0.04 (rank 51 % within c5).
- **rs2395029** (HLA-B*57:01 tag; in_NMF_set = True): best component **c1** (H = 0.47, rank 32 % within c1); H_c5 = 0 (rank 64 % within c5; not in the top 5 % of c5).
- **rs2244546** (Lange 2013 tag; in_NMF_set = True): best component **c3** (H = 0.10, rank 59 % within c3), with smaller secondary loadings on c0 and c6; H_c5 = 0 (not in the top 5 % of c5).
- **rs11509487** (MICB 9-bp indel, paper-1 anchor; filtered as non-biallelic; nearest in-NMF SNV at 31,475,127, 43 bp away): nearest SNV's argmax is c2 (H = 0.49, **rank 12 % within c2**); H_c5 = 0.37 (rank 34 % within c5; not in the top 5 % of c5).

None of the four manuscript-relevant tag SNPs (or their nearest in-NMF SNVs) fall in the **top 5 %** of any component's H loading distribution, and in particular none of them fall in the top 5 % of **c5** (the HLA-B / HLA-C-variable axis). This supports the manuscript claim that the c5 axis is not captured by previously published HCV-HCC tag SNPs in the MICA / HCP5 / MICB region — c5 represents a haplotypic background that is uncharacterised in the existing tag-SNP literature.

**Table S9 — LIRI-JP regression sensitivities and multiple-testing.**
Source: `supp_6_liri_sensitivity.py`. Two CSV files:

- `table_S9_liri_sensitivity.csv` — per-row analysis label (Primary n=117, S1 whitespace-trim n=120, S2 + HBV+HCV→HCV n=122, S3 c5 × viral interaction terms), c5 β, SE, 95 % CI, p, R². With the per-analytic-sample z-standardisation harmonised across scripts (ddof = 0 within the filtered sample), Primary β = −0.193 (p = 5.3 × 10⁻⁴) matches Table S1 section A row M3 exactly. S1 (n=120) β = −0.185 and S2 (n=122) β = −0.187 confirm the c5 → HLA-C ↓ effect is robust to whitespace-trim recovery and HBV+HCV co-infection re-assignment. Interaction terms NS (HCV p = 0.24, NBNC p = 0.78).
- `table_S9b_multiple_testing.csv` — Table S1 section C (7 gene outcomes under M3 full model) augmented with BH-FDR and Bonferroni adjusted p-values. HLA-C survives both corrections (adjusted p = 0.0037); the other six outcomes do not.

---

**Table S10 — Cross-tissue NES direction analysis (Axis I + Axis II).**
Source: `axis_sign_coherence_matrix.py`. Companion to Figure 3 panel (c).

- `table_S10_axis_sign_coherence.csv` — for every (component × tissue × gene) combination in {c4, c5, c6} × {Whole_Blood, Liver, Cells_EBV-transformed_lymphocytes (LCL), Lung, Skin_Sun_Exposed_Lower_leg, Colon_Transverse} × {MICA, HLA-B, HLA-C}, reports n_eqtl_in_window, n_top_overlap, n_pos, n_neg, sign_coherence, mean_NES, median_NES. Axis I cells (c4 / c6 × MICA): every one of 12 reported cells shows sign coherence = 1.00 in the negative direction (mean NES range −0.35 in Whole Blood to −0.72 in Colon Transverse; combined 2,731 SNV-tissue observations, n_pos = 0). Axis II cells (c5 × HLA-B, c5 × HLA-C): all 4–5 observed tissues per gene show sign coherence = 1.00 (HLA-B positive, HLA-C negative). Off-axis cells (e.g. c4 × HLA-C, c6 × HLA-C, c5 × MICA) included for completeness and show appreciably lower sign coherence (0.54–0.97), demonstrating that NMF cleanly separates the two axes by primary gene target rather than capturing a single hand-of-MHC haplotype.

---

## LIRI-JP analytic sample-size note (n = 122 vs n = 117)

LIRI-JP included **122** donors with both germline genotype (rs2596542 carriers + c5 top-signature SNVs available in `LIRI-JP_MHC_germline.vcf.gz`) and tumour RNA-seq (EGA-derived TPM matrix; gene symbols indexed). Univariate analyses (raw c5 score vs HLA-C, etc.) and the c5 tertile boxplot in main Figure 5 (a) use n = 122.

All **multivariable regression models** (Table S1 sections A–C, Supplementary Figures S4 and S5) require complete data on c5 score, expression, T_stage, viral_status (∈ {HBV, HCV, NBNC}), age, sex, and HLA-A / immune-marker TPM. Five donors are dropped from the complete-case subset because of missing T_stage or non-standard viral_status entries (`HBV, HCV` co-infection, n = 2; whitespace-suffix entries `HBV `, `HCV `, n = 3). The complete-case n is therefore **117**.

Manuscript wording (insert into Methods):

> *LIRI-JP included 122 donors with matched genotype and RNA-seq data; complete-case multivariable regression models including all covariates (T_stage, viral_status, age, sex, immune-marker expression, log₂(HLA-A)) were performed in 117 donors after exclusion of 5 donors with missing T_stage or non-standard viral_status entries.*

---

## Excluded from supplementary (decisions log)

- **Survival analysis** (`liri_c5_survival.py`, `results/liri_c5_survival_cox.csv`, `results/figure_liri_c5_survival.{pdf,png}`) — DROPPED on 2026-05-08. n = 117 with only 16 disease-specific death events, too underpowered for primary AJHG submission. Files kept in repo for reproducibility and for re-analysis if a replication cohort becomes available, but not cited, summarised, or supplied as Supplementary in the manuscript.
- **Internal PASS/FAIL summary** (`results/supp_combined_summary.csv`) — internal QA artifact, not for publication. The thresholds it encodes (matched corr > 0.7, cophenetic-within-0.05, score corr > 0.8) are post-hoc rules-of-thumb chosen during analysis design and should not be presented as pre-registered tests.

---

## File-system layout

```
supplementary/
├── SUPPLEMENTARY_INDEX.md
├── figures/
│   ├── figure_S1_lda_vs_nmf.{pdf,png}
│   ├── figure_S2_cophenetic_k.{pdf,png}
│   ├── figure_S3_seed_stability.{pdf,png}
│   ├── figure_S4_liri_weighted.{pdf,png}
│   ├── figure_S5_liri_hlaa_diagnostic.{pdf,png}
│   ├── figure_S6_liri_sensitivity.{pdf,png}
│   ├── figure_S7_comt_flipflop.{pdf,png}
│   └── _archive/                                     ← legacy figures (not cited from manuscript)
│       └── figure_2_branch_composition.{pdf,png}     ← old k = 2..5 hierarchical-clustering Figure 2
└── tables/
    ├── table_S1_regression_models.csv
    ├── table_S1d_VIF.csv
    ├── table_S2_component_top_snvs.csv
    ├── table_S3_matched_permutation.csv
    ├── table_S4_gtex_tissue_eqtl_NES.csv
    ├── table_S5_lda_vs_nmf.csv
    ├── table_S6_nmf_rank_cophenetic.csv
    ├── table_S7_nmf_seed_stability.csv
    ├── table_S7b_nmf_seed_stability_per_seed.csv
    ├── table_S8_tag_snp_component_mapping.csv
    ├── table_S8_tag_snp_component_mapping_wide.csv
    ├── table_S9_liri_sensitivity.csv
    ├── table_S9b_multiple_testing.csv
    └── table_S10_axis_sign_coherence.csv
```
