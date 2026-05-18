# Figure and table inventory — final (pre-AJHG submission)

Date: 2026-05-15
Repo: `~/paper1-proxy-topology/mica_flipflop/`
Supersedes: `docs/FIGURE_NUMBERING.md` Figure-2 row (which still describes the legacy k=2..5 branch-composition figure; see `reports/addon_audit_final_figures_tables.md` issue F2-1).

Five main figures + seven supplementary figures + ten supplementary tables.

---

## Main figures

| # | File (PDF/PNG) | Source script | One-line narrative |
|---|---|---|---|
| **Fig 1** | `results/figure_1_lange_replication.{pdf,png}` | `figure_1_lange_replication.py` | Cross-population signed-LD reversal — Lange 2013 Fig. 1A reproduced in 1000G phase 3 across 10 sub-pops (5 EAS + 5 EUR), harmonised to Lange-coded alleles for rs2596542 × rs2244546 and rs2596542 × rs9275572. |
| **Fig 2** | `figures/main/figure_2_nmf_k8_population_composition.{pdf,png}` | `figure_2_nmf_k8_population_composition.py` | k = 8 NMF component composition among the 2,111 rs2596542-T carrier haplotypes, across the 5 EAS + 5 EUR Lange-1A sub-populations. Panel (a) mean component loading; panel (b) row-normalised composition. **c4 + c6 = Axis I (MICA, EAS-enriched)**; **c5 = Axis II (HLA-B ↑ / HLA-C ↓, EUR-enriched)**. |
| **Fig 3** | `results/figure_3_component_eqtl_axes.{pdf,png}` | `eqtl_axis_test.py` (panels A, B), `axis_sign_coherence_matrix.py` (panel C) | Component × eQTL regulatory-axis annotation. Each NMF component maps to a distinct eQTL regulation axis (GTEx Whole_Blood). c2 / c4 / c6 = MICA-stable axes; c5 = HLA-B-variable axis (HLA-B↑, HLA-C↓, HCG27↑, MICB↑). Panel (c) cross-tissue NES sign-coherence matrix references `Supplementary Table S10`. |
| **Fig 4** | `results/figure_4_robustness.{pdf,png}` | `permutation_all_components.py` | Matched-permutation robustness. Per (NMF component × test set) MAF×distance matched permutation (10,000 perms) confirms c5 enrichment is not driven by SNV density or anchor proximity. |
| **Fig 5** | `results/figure_5_liri_panel.{pdf,png}` | `figure_5_liri_panel.py` | LIRI-JP HCC validation (n = 117 with full covariates). (a) c5 score tertile vs log₂(HLA-C TPM+1); (b) forest plot of c5 → HLA-C ↓ under progressive adjustment M0 → M3. β = −0.193, p = 5.3 × 10⁻⁴ in M3. |

---

## Supplementary figures

| # | File (PDF/PNG) | Source script | One-line narrative |
|---|---|---|---|
| **Fig S1** | `supplementary/figures/figure_S1_lda_vs_nmf.{pdf,png}` | `supp_1_lda_comparison.py` | NMF vs LDA on the same haplotype matrix. Hungarian-matched component correlations c4 / c5 / c6 = 0.99 / 0.87 / 0.96. LDA preserves MICA depletion and C-flip enrichment but attenuates the localised HLA-B *cis*-eQTL signal (NMF fold 2.14 → LDA fold 0.66). |
| **Fig S2** | `supplementary/figures/figure_S2_cophenetic_k.{pdf,png}` | `supp_2_cophenetic_k.py` | NMF rank selection by consensus cophenetic correlation across k = 2..12, 20 random seeds per rank. cophenetic max 0.998 at k = 2; k = 8 cophenetic 0.963 (within 0.05 of max). |
| **Fig S3** | `supplementary/figures/figure_S3_seed_stability.{pdf,png}` | `supp_3_seed_stability.py` | NMF seed reproducibility at k = 8 across 50 random initialisations. Matched-corr means c4 / c5 / c6 = 0.98 / 0.96 / 0.99; MICA-enriched label recovered in 94 % of seeds, HLA-B-enriched in 100 %. |
| **Fig S4** | `supplementary/figures/figure_S4_liri_weighted.{pdf,png}` | `supp_4_liri_weighted.py` | LIRI-JP c5 score weighting sensitivity. Unweighted vs H-loading-weighted vs NNLS-projected definitions in n = 117 full-covariate analytic sample. All three retain c5 → HLA-C ↓ (β = −0.193 / −0.160 / −0.175; all P < 0.005). |
| **Fig S5** | `supplementary/figures/figure_S5_liri_hlaa_diagnostic.{pdf,png}` | `supp_5_liri_hlaa_conditioning.py` | HLA-A adjustment diagnostic. c5 ⊥ HLA-A (r = −0.04, p = 0.64); HLA-A captures 70 % of HLA-C variance; SE shrinkage M2 → M3 = −30 %. HLA-A is a precision-improving covariate, not a confounder. |
| **Fig S6** | `supplementary/figures/figure_S6_liri_sensitivity.{pdf,png}` | `supp_6_liri_sensitivity.py` | LIRI-JP HLA-C regression sensitivity (analytic-sample S1 / S2 / S3 alternatives) and multiple-testing context (BH-FDR + Bonferroni across the 7-gene outcome family). Only HLA-C survives both corrections (adjusted p = 0.0037). |
| **Fig S7** | `supplementary/figures/figure_S7_comt_flipflop.{pdf,png}` | `figure_S7_comt_flipflop.py` | COMT Val158Met (rs4680) reference-panel replication and NMF extension. (A) Population-conditional signed LD and asymmetry C for rs4680 × rs2097603; (B) k = 4 NMF mixture conditional on rs4680 (V/M); (C) population-level COMT NMF mixture across 5,008 haplotypes. Proof-of-concept that carrier-set decomposition generalises beyond the MHC anchor. Not cited from main Results. |

**Archived (unnumbered, not cited):** `supplementary/figures/_archive/figure_2_branch_composition.{pdf,png}` and `supplementary/figures/_archive/figure_S_coarse_carrier_space_k2_k5.{pdf,png}` — legacy main Figure 2 (k = 2..5 hierarchical-clustering branch composition); kept for reproducibility but superseded by the k = 8 NMF view (Fig 2 above).

---

## Supplementary tables

| # | File | Source script | One-line narrative |
|---|---|---|---|
| **S1** | `supplementary/tables/table_S1_regression_models.csv` (+ `table_S1d_VIF.csv`) | `make_supplementary_package.py` | LIRI-JP regression and sensitivity analyses. Three sections: (A) HLA-C ~ c5 progressive adjustment M0–M3; (B) HLA-C M3 with three c5 score variants; (C) M3 full model for each of 7 gene outcomes. Companion S1d VIF table. |
| **S2** | `supplementary/tables/table_S2_component_top_snvs.csv` | `make_supplementary_package.py` | Top-signature SNVs defining NMF components (8 × 355 = 2,840 rows). Columns: component, rank_within_component, position_grch37, H_loading (4 d.p.), gene_region, C_flip, class_flip, eQTL_{MICA, MICB, HLA-B, HLA-C, HCG27, HCP5}. Anchor rs2596542 excluded by construction. |
| **S3** | `supplementary/tables/table_S3_matched_permutation.csv` | `permutation_all_components.py` | MAF × distance matched permutation results (10,000 perms) per (component × test set). Test sets: r-flip / C-flip / class-flip / MICA / MICB / HLA-B / HLA-C eQTL. |
| **S4** | `supplementary/tables/table_S4_gtex_tissue_eqtl_NES.csv` | `eqtl_axis_test.py`, `liver_replication.py`, `external_tissue_replication.py` | GTEx tissue eQTL NES summary per (component × tissue × target gene). Two analysis blocks: 8-component primary (Whole_Blood + Liver) and c5 cross-tissue replication (LCL + Lung + Skin). |
| **S5** | `supplementary/tables/table_S5_lda_vs_nmf.csv` | `supp_1_lda_comparison.py` | LDA vs NMF c5-analogue fold enrichment on 6 test sets (Hungarian-matched). |
| **S6** | `supplementary/tables/table_S6_nmf_rank_cophenetic.csv` | `supp_2_cophenetic_k.py` | NMF rank cophenetic and reconstruction-error table for k = 2..12 (20 seeds per rank). |
| **S7** | `supplementary/tables/table_S7_nmf_seed_stability.csv` (+ `table_S7b_nmf_seed_stability_per_seed.csv`) | `supp_3_seed_stability.py` | NMF seed stability per component (50 random seeds). |
| **S8** | `supplementary/tables/table_S8_tag_snp_component_mapping.csv` (+ `table_S8_tag_snp_component_mapping_wide.csv`) | `tag_snp_component_map.py`, `make_supplementary_package.py` | Tag SNP → NMF component mapping. 4 manuscript-relevant tag SNPs × 8 components. None of {rs2596542 (anchor), rs2244546, rs2395029, rs11509487} fall in the top 5 % of c5 — the c5 axis is uncharacterised by existing tag-SNP literature. |
| **S9** | `supplementary/tables/table_S9_liri_sensitivity.csv` (+ `table_S9b_multiple_testing.csv`) | `supp_6_liri_sensitivity.py` | LIRI-JP sample-definition sensitivity (Primary 117 / S1 120 / S2 122) and 7-outcome BH-FDR + Bonferroni multiple-testing correction. Only HLA-C survives correction. |
| **S10** | `supplementary/tables/table_S10_axis_sign_coherence.csv` | `axis_sign_coherence_matrix.py` | Axis sign-coherence summary (Axis I + Axis II) across 6 GTEx tissues × {MICA, HLA-B, HLA-C}. Axis I cells (c4/c6 × MICA): sign coherence = 1.00 (mean NES < 0). Axis II cells (c5 × HLA-B / HLA-C): sign coherence = 1.00 in expected directions. |

---

## Companion working tables (not supplementary)

| File | Purpose |
|---|---|
| `tables/figure_2_component_population_summary.csv` | Long-format companion to revised Figure 2 (26 pops × 8 components × `mean_loading`, `carrier_fraction`, `frac_nonzero`, `n_haps`, `component_role`). Source: `figure_2_nmf_k8_population_composition.py`. |

---

## Conventions

- **Coordinates**: GRCh37 (1000G phase 3 native) in main figures. dbSNP / GTEx (GRCh38) requires offset chr6 +32,223 bp — handled in `liver_replication.py` and `external_tissue_replication.py`.
- **Allele coding**: dbSNP forward-strand ALT = 1 throughout. Figure 1 only is harmonised to Lange 2013 coded alleles (rs2596542-T, rs2244546-G, rs9275572-A) per `docs/rs2596542_allele_strand_build_mapping.md` §4b.
- **Anchor exclusion**: rs2596542 is removed from the NMF input matrix after carrier selection (zero variance among the rs2596542-T carriers).
- **Top-signature threshold**: top 5 % of H[c, ·] per component = 355 SNVs (Table S2). Sensitivity to 1 / 2.5 / 10 / 15 % thresholds reported in Supplementary Methods.
- **LIRI-JP analytic samples**: 122 with matched genotype + RNA-seq; 117 complete-case for multivariable models (5 dropped for missing T_stage or non-standard viral_status entries).
- **c5 score**: unweighted "c5 signature burden" score (count of available top-signature SNVs in LIRI-JP genotype) used as primary; H-loading-weighted and NNLS-projected scores reported as sensitivity analyses (Table S1 §B, Fig S4).

---

## Outstanding decisions — resolved 2026-05-15

- **U-2** (legacy k = 2..5 figure numbering) — resolved: not assigned a Supplementary Figure number. Archived at `supplementary/figures/_archive/`.
- **U-3** (supplementary figure dir consolidation) — resolved: `supplementary/figures/` is the canonical directory for all supplementary figure files. The 26-pop Fig 2 companion remains alongside the main Fig 2 render under `figures/supplementary/figure_2_nmf_k8_26pop.{pdf,png}` because it is a direct companion produced by the Figure 2 source script.
- **U-4** (`tables/figure_2_component_population_summary.csv`) — resolved: treated as Figure 2 source data, not a numbered Supplementary Table.

## Still outstanding

- **U-1**: Where is the manuscript DOCX / master source? Not yet stored in this repo.

---

End of inventory.
