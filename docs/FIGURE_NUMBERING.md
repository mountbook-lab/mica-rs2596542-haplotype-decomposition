# Manuscript figure numbering — locked

Last updated: 2026-05-08

Five main figures + supplementary panels. Each row links the manuscript
figure number to (i) the source script, (ii) the generated PDF/PNG with
the manuscript-numbered filename, (iii) the underlying CSV / parquet
data, and (iv) the one-line narrative.

## Main figures

| # | File (PDF/PNG) | Source script | Narrative |
|---|---|---|---|
| **Fig 1** | `results/figure_1_lange_replication.{pdf,png}` | `figure_1_lange_replication.py` | Lange 2013 Fig. 1A reproduced in 1000G phase 3 across 10 sub-pops, harmonized to Lange-coded alleles. Population-conditional sign of LD between rs2596542 (chr6:31,366,595 GRCh37, T) and two HCV-HCC tag candidates (rs2244546-G chr6:31,435,833; rs9275572-A chr6:32,678,999). |
| **Fig 2** | `results/figure_2_branch_composition.{pdf,png}` | `cluster_anchor_haplotypes.py` (renders `figure_branch_composition`) | Anchor-carrier haplotype branch composition across 1000G sub-pops (rs2596542-T carriers; k=2..5 hierarchical clustering on Hamming distance, ±250 kb window). Mechanistic explanation of the Lange flip: H1 frequency varies ~7-fold between JPT (0.55) and FIN (0.08) **among carriers of the same anchor allele**, so the same tag SNP captures different functional backgrounds in different populations → opposite-signed correlations with nearby SNVs. |
| **Fig 3** | `results/figure_3_component_eqtl_axes.{pdf,png}` | `eqtl_axis_test.py` (renders `figure_component_eqtl_axis`) | Each NMF component maps to a distinct eQTL regulation axis (GTEx Whole_Blood). c2 / c4 / c6 = MICA-stable axes; c5 = HLA-B-variable axis (HLA-B↑, HLA-C↓, HCG27↑, MICB↑). |
| **Fig 4** | `results/figure_4_robustness.{pdf,png}` | `permutation_all_components.py` (renders `figure_permutation_all_components`) | Per-component MAF×distance matched permutation (10,000 permutations) confirms c5 enrichment is not driven by SNV density or anchor-proximity. Heatmaps: z-score (Panel A) and fold vs matched (Panel B) for r-flip / C-flip / class-flip / HLA-B / HLA-C / MICB / MICA test sets. |
| **Fig 5** | `results/figure_5_liri_panel.{pdf,png}` | `figure_5_liri_panel.py` | LIRI-JP HCC validation (n=117 with full covariates). (a) c5 score tertile vs log2(HLA-C TPM+1), etiology-colored — Jonckheere trend p=0.008. (b) Forest plot: c5 → HLA-C ↓ robust under T-stage / viral / age / sex / immune-PC1 / HLA-A adjustment. |

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
| `figure_liri_c5_survival.{pdf,png}` | LIRI Cox + KM survival by c5 tertile |
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
