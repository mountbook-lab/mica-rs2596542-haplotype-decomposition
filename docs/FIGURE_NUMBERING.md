# Manuscript figure numbering — locked

Last updated: 2026-05-15

Five main figures + supplementary panels. Each row links the manuscript
figure number to (i) the source script, (ii) the generated PDF/PNG with
the manuscript-numbered filename, (iii) the underlying CSV / parquet
data, and (iv) the one-line narrative.

## Main figures

| # | File (PDF/PNG) | Source script | Narrative |
|---|---|---|---|
| **Fig 1** | `results/figure_1_lange_replication.{pdf,png}` | `figure_1_lange_replication.py` | Lange 2013 Fig. 1A reproduced in 1000G phase 3 across 10 sub-pops, harmonized to Lange-coded alleles. Population-conditional sign of LD between rs2596542 (chr6:31,366,595 GRCh37, T) and two HCV-HCC tag candidates (rs2244546-G chr6:31,435,833; rs9275572-A chr6:32,678,999). Short black horizontal ticks at the JPT and CEU bars mark the values reported by Lange et al. (2013); minor deviations from this study's 1000G phase 3 estimates (e.g. for rs2244546 JPT +0.49 vs Lange +0.59; for rs9275572 JPT +0.48 vs Lange +0.27) reflect differences in sample composition between Lange's HapMap3-derived cohort and 1000G phase 3 JPT, and (per Suppl-Fig 1) the resulting allele-frequency shifts that modulate \|r\| while preserving the sign. The rs2244546 CEU bar coincides numerically with Lange's reported CEU value (both −0.19), so the reference tick is visually fused with the top edge of the bar. |
| **Fig 2** | `figures/main/figure_2_nmf_k8_population_composition.{pdf,png}` | `figure_2_nmf_k8_population_composition.py` (consumes `results/nmf_W_26_k8.parquet` + `results/_supp_hap_meta.parquet` written by `decompose_branches_26.py`) | **Figure 2. k = 8 NMF component composition among rs2596542-T carrier haplotypes in EAS and EUR populations.** Non-negative matrix factorization (k = 8) on the 2,111 rs2596542-T carrier haplotypes × 7,116 SNV matrix, shown across the 5 EAS + 5 EUR sub-populations of Lange 2013 Fig. 1A. Panel (a) per-population mean NMF loading; panel (b) row-normalized component composition. **c4 + c6 = Axis I (MICA, EAS-enriched)**; **c5 = Axis II (HLA-B↑ / HLA-C↓, EUR-enriched)**. A 26-population companion is provided as Supplementary Figure 2 (`figures/supplementary/figure_2_nmf_k8_26pop.{pdf,png}`); the legacy k = 2..5 hierarchical-clustering branch figure is archived at `supplementary/figures/_archive/figure_2_branch_composition.{pdf,png}` and is no longer cited in the main text. |
| **Fig 3** | `results/figure_3_component_eqtl_axes.{pdf,png}` | `eqtl_axis_test.py` (renders `figure_component_eqtl_axis`) | Each NMF component maps to a distinct eQTL regulation axis (GTEx Whole_Blood). c2 / c4 / c6 = MICA-stable axes; c5 = HLA-B-variable axis (HLA-B↑, HLA-C↓, HCG27↑, MICB↑). Panel A: fold enrichment of each component's top-5 % SNV set against per-gene eQTL targets (color scale capped at 4). HCP5 fold-enrichment values are not displayed in panel A because the expected background overlap exceeds ~95 % within the ±200 kb window (HCP5 has dense cis-eQTL coverage across the MHC), rendering fold enrichment uninformative. Panel B: hypergeometric significance, plotted as −log₁₀ p with the color scale capped at 15 (the printed cell value is the true −log₁₀ p; deep-blue cells are saturated, not capped numerically). Cells labelled `>308` are exact-test p-values that fall below the float64 underflow threshold (−log₁₀ p > 308, i.e. p < 10⁻³⁰⁸) under the hypergeometric independence assumption — this independence assumption breaks down in HCP5's high-density regime, which is why panel A suppresses the corresponding fold values. |
| **Fig 4** | `results/figure_4_robustness.{pdf,png}` | `permutation_all_components.py` (renders `figure_permutation_all_components`) | Per-component MAF×distance matched permutation (10,000 permutations) confirms c5 enrichment is not driven by SNV density or anchor-proximity. Heatmaps: z-score (Panel A) and fold vs matched (Panel B) for r-flip / C-flip / class-flip / HLA-B / HLA-C / MICB / MICA test sets. |
| **Fig 5** | `results/figure_5_liri_panel.{pdf,png}` | `figure_5_liri_panel.py` | LIRI-JP HCC validation (n=117 with full covariates). (a) c5 score tertile vs log2(HLA-C TPM+1), etiology-colored. Trend across tertiles: Jonckheere–Terpstra two-sided p = 0.008; low vs high Mann–Whitney U p = 0.008; continuous c5 vs HLA-C Spearman ρ = −0.28, p = 0.002 (statistics moved out of the panel to avoid overlap with the low-tertile scatter). (b) Forest plot of the c5 coefficient on log2(HLA-C) under progressively stricter adjustment, using the M0–M3 numbering of Table S1 section A (M0 = c5 only, omitted from the forest plot as it is not an adjustment step): M1 = + clinical (T-stage + viral + age + sex); M2 = + immune-PC1; M3 = + log2(HLA-A) [strict]. The independent "+ PTPRC" step shown in earlier drafts was removed because PTPRC is one of the six markers defining immune-PC1, so it would double-count rather than form an independent stage. c5 → HLA-C ↓ robust through M3 (β = −0.193, p = 5.3 × 10⁻⁴ in M3; matches Table S1 section A row M3). |
| **Fig 3 update (panel c)** | `results/figure_component_eqtl_axis.{pdf,png}` (panel c) | `axis_sign_coherence_matrix.py` (builds the matrix), rendered by `eqtl_axis_test.py` | Cross-tissue direction analysis of the two NMF axes against gene-level cis-eQTL NES (GTEx v8). Axis I rows (c4 × MICA, c6 × MICA): all 12 panel cells show sign coherence = 1.00 with mean NES < 0 (1,010 + 1,721 = 2,731 SNV-tissue observations, all NES < 0), confirming unidirectional MICA-down regulation across six tissues — Whole Blood, Liver, EBV-lymphocytes (LCL), Lung, Skin Sun-Exposed, Colon Transverse. Axis II row 1 (c5 × HLA-B): all 4 tissues with overlap show sign coherence = 1.00 with mean NES > 0 (HLA-B-up). Axis II row 2 (c5 × HLA-C): all 5 tissues with overlap show sign coherence = 1.00 with mean NES < 0 (HLA-C-down). LCL has no overlapping HLA-B eQTL records ("n/a"). The accompanying Supplementary Table S10 (`table_S10_axis_sign_coherence.csv`) reports the full long-format breakdown including the off-axis cells (c4 × HLA-C, c6 × HLA-C, c4 × HLA-B, c6 × HLA-B, c5 × MICA) that exhibit lower sign coherence (range 0.54–0.99), supporting the conclusion that NMF cleanly separates the two axes by primary gene target. |

## Supplementary panels

| # | File (PDF/PNG) | Source script | Narrative |
|---|---|---|---|
| Suppl-Fig 1 | `results/figure_1_supplementary_AF.{pdf,png}` | `figure_1_supplementary_AF.py` | Allele-frequency context for Lange Fig. 1A magnitude differences. Anchor (rs2596542-T) p across 10 sub-pops with Lange's reported HapMap3 values (0.329 JPT, 0.284 CEU); partner (rs2244546-G, rs9275572-A) p across 10 sub-pops. Magnitude discrepancies between HapMap3 and 1000G are AF-driven, not coding errors. |

## Other figures (candidate supplementary; not yet renumbered)

These were generated during the analysis arc and remain in
`results/figure_*.{pdf,png}` for reference. If any are promoted to a
supplementary slot, copy with the `figure_S{N}_` prefix.

| Source filename | Topic |
|---|---|
| `figure_branch_pca.{pdf,png}` | PCA on anchor-T carriers — visual sanity for branch existence |
| `figure_branch_composition.{pdf,png}` | source for Fig 2 (renamed copy) |
| `figure_umap_leiden_26.{pdf,png}` | UMAP+Leiden clustering on 26-pop anchor-T carriers |
| `figure_nmf_components_26.{pdf,png}` | NMF k=8 SNV signature matrix on 26-pop carriers — supplementary technical detail of which SNVs make up each component |
| `figure_nmf_components.{pdf,png}` | k=8 NMF on 10-pop subset (precursor) |
| `figure_component_annotation.{pdf,png}` | gene-region / tag-SNP annotation overlay per component |
| `figure_tag_snp_component.{pdf,png}` | rs2244546 / rs2395029 / rs11509487 → component map |
| `figure_eqtl_direction.{pdf,png}` | sign coherence of NES per (component, target gene) |
| `figure_external_tissue_replication.{pdf,png}` | LCL / Lung / Skin GTEx replication of c5 |
| `figure_liver_replication.{pdf,png}` | GTEx Liver replication of c5 |
| `figure_k_sensitivity.{pdf,png}` | k = 4..12 robustness of c5 enrichment |
| `figure_threshold_sensitivity.{pdf,png}` | top-1% / 5% / 10% threshold sensitivity |
| `figure_ld_metric_comparison.{pdf,png}` | r² / signed_r / D / D' / Δ / C / M / abs_r flip-detection comparison |
| `figure_simulation_ground_truth.{pdf,png}` | NMF correctly recovers known branches in synthetic data |
| `figure_c5_permutation.{pdf,png}` | c5-only matched permutation null (precursor to Fig 4) |
| `figure_MC_plane.{pdf,png}` | M/C plane for the 2-SNP example |
| `figure_region_MC.{pdf,png}` | region-wide M/C scan |
| `figure_signed_r_vs_topology.{pdf,png}` | signed r × topology class (raw dbSNP coding — superseded by Fig 1) |
| `figure_lange_replication_grch37.{pdf,png}` | preliminary Lange replication (superseded by figure_1_lange_replication) |
| `figure_flip_type_comparison.{pdf,png}` | r-flip / C-flip / class-flip cross-tab |
| `figure_within_branch_scan.{pdf,png}` | within-branch carrier-set topology scan |
| `figure_afr_branches.{pdf,png}` | AFR-private branch deepdive (extended analysis) |
| `figure_liri_c5_validation.{pdf,png}` | LIRI 7-gene scatter (precursor to Fig 5) |
| `figure_liri_c5_clinical.{pdf,png}` | LIRI clinical correlations (c5 score / HLA-C × age, gender, viral, T-stage) |
| `figure_liri_c5_survival.{pdf,png}` | LIRI Cox + KM survival by c5 tertile — **EXCLUDED from AJHG submission package** (n = 16 events; underpowered). Kept in `results/` for reproducibility only; not cited in main text, manuscript-revisions doc, or `supplementary/`. See `supplementary/SUPPLEMENTARY_INDEX.md` "Excluded from supplementary" section. |
| `figure_liri_c5_tertile.{pdf,png}` | original (b)/(d) labels — superseded by figure_5_liri_panel |
| `figure_component_eqtl_axis.{pdf,png}` | source for Fig 3 (renamed copy) |
| `figure_permutation_all_components.{pdf,png}` | source for Fig 4 (renamed copy) |

## Convention notes

- **Coordinates**: all positions in main figures use **GRCh37** (1000G
  phase 3 native coordinate system). When citing dbSNP / GTEx (which
  use GRCh38), the offset chr6 +32,223 bp must be applied — handled
  automatically in `liver_replication.py` and `external_tissue_replication.py`.
- **Allele coding**: default coding throughout the code is dbSNP
  forward-strand ALT = 1. For Fig 1 only, signed correlations are
  harmonized to Lange 2013 coded alleles (rs2596542-T, rs2244546-G,
  rs9275572-A) per `rs2596542_allele_strand_build_mapping.md` §4b.
- **Renumbered files are copies, not moves**. Original
  `figure_nmf_components_26.{pdf,png}` etc. are untouched — analysts
  who run the upstream scripts will continue to see the original
  filenames in `results/`.
