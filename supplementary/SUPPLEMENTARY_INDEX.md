# Supplementary materials index — pre-AJHG submission

Date: 2026-05-15.  Project: `mica-rs2596542-haplotype-decomposition`.

All numbers in the captions below are taken directly from the corresponding table CSV. Source scripts are in the parent repo.

---

## Supplementary Figures (4)

**Figure S1 — Comparison of NMF and LDA decomposition on the same rs2596542-T carrier haplotype matrix.**
Source: `supp_1_lda_comparison.py`. k = 8 LatentDirichletAllocation fitted on the same 2,111 × 7,116 haplotype matrix; Hungarian matching of LDA components to NMF components via H-row Pearson correlation. Panel (a) — matched correlation per NMF component; the dashed red line marks the matching threshold r = 0.7. Panel (b) — LDA-c5 analogue eQTL / C-flip enrichment alongside NMF-c5 (c5-like enrichment pattern is *partially retained* under LDA); cells annotated `n/a` indicate zero overlap (fold = 0) between the matched component's top-5 % SNV set and the gene's eQTL set in the universe (no MICA eQTL overlap in either method's c5-like component). **Key values**: matched corr c4 / c5 / c6 = 0.99 / 0.87 / 0.96; LDA-c5 preserves MICA depletion (fold = 0) and C-flip enrichment (fold = 1.39) but attenuates the localized HLA-B *cis*-eQTL signal (NMF fold 2.14 → LDA fold 0.66). LDA recovered the major component structure but smoothed the localized HLA-B regulatory signature, supporting NMF as the appropriate method for sparse, localized component-defining SNV signatures.

**Figure S2 — NMF rank selection by consensus cophenetic correlation.**
Source: `supp_2_cophenetic_k.py`. k = 2..12 with 20 random initializations per rank; consensus matrix C[i,j] = (fraction of seeds where haplotypes i and j share dominant component); cophenetic correlation = Pearson r between cophenetic distances on the average-linkage dendrogram and (1 − C). Panel (a) — cophenetic correlation and Brunet-2004 dispersion vs k; (b) — mean reconstruction error ‖X − WH‖_F. **Key values**: cophenetic max = 0.998 at k = 2; k = 8 cophenetic = 0.963 (within 0.05 of max).

**Figure S3 — NMF seed reproducibility at k = 8.**
Source: `supp_3_seed_stability.py`. 50 independent NMF runs at k = 8 (init = random). Hungarian-matched H-row Pearson correlations (anchored to seed 0); recovery of MICA-enriched (fold > 3) and HLA-B-enriched (fold > 1.5) component labels per seed. Panel (a) — boxplot of matched corr per component; (b) — MICA / HLA-B fold across 49 matched seeds. **Key values**: c4 / c5 / c6 mean matched corr = 0.98 / 0.96 / 0.99; MICA-enriched recovered in 46/49 (94 %) seeds; HLA-B-enriched recovered in 49/49 (100 %) seeds.

**Figure S4 — Reference-panel replication and NMF extension of the COMT Val158Met flip-flop example.**
Source: `figure_S4_comt_flipflop.py` (consumes `data/comt_flipflop_nmf/tables/comt_pair_metrics_by_pop.tsv`, `comt_component_mixture_by_rs4680_status.tsv`, `comt_component_mixture_by_pop.tsv`). (A) Population-conditional signed LD and asymmetry component C for rs4680 (Val158Met) and rs2097603 across 1000 Genomes populations. Blue and brown indicate positive and negative signed values, respectively. (B) Row-normalized k = 4 NMF component mixture conditional on rs4680 allele status, with V indicating Val/REF carriers and M indicating Met/ALT carriers. (C) Population-level COMT NMF component mixture across all 5,008 haplotypes. Compared with the MHC rs2596542 locus, the COMT decomposition is simpler and more population-like, supporting broader applicability of carrier-set decomposition while illustrating locus-dependent differences in functional interpretability. COMT is included as an external reference-panel proof-of-concept and is not cited in the main Results.

## Supplementary Tables (9)

**Table S1 — Top-signature SNVs defining NMF components in rs2596542-T carrier haplotypes.**
Built by `make_supplementary_package.py` from `results/nmf_H_26_k8.parquet` and `results/region_per_partner_26.csv`. Top-signature SNVs were defined as variants with H-matrix loadings in the top 5 % of each component's H[c, ·] distribution (i.e., at or above the 95th percentile), yielding 355 SNVs per component (8 × 355 = 2,840 rows). The table lists the component-defining SNV sets used for LD-reversal enrichment, *cis*-eQTL annotation, MAF×distance matched-permutation analysis (Supplementary Table S2), and tag-SNP component mapping (Supplementary Table S7). The anchor SNV rs2596542 was removed from the NMF input matrix after carrier selection because all retained haplotypes carry the rs2596542-T allele at this site, leaving it with zero variance in the analytic matrix. Columns: `component` (c0..c7); `rank_within_component` (1 = highest H_loading, descending); `position_grch37` (chr6 1000G phase 3 coordinate); `H_loading` (rounded to 4 decimal places); `gene_region` (MICA / MICB / HLA-B / HLA-C / HCG27 / HCP5 / 3'_of_MICB / 5'_of_HCG27 / intergenic_class_I); `C_flip`, `class_flip` (boolean LD-reversal flags); `eQTL_<gene>` (boolean indicating whether the SNV is a GTEx Whole_Blood *cis*-eQTL for that gene). **Rows**: 2,840 (8 components × 355).

**Table S2 — Matched permutation for component-level enrichment.**
Source: `permutation_all_components.py`. Per (NMF component × test set) MAF×distance matched permutation, 10,000 perms. Test sets: r-flip / C-flip / class-flip / MICA / MICB / HLA-B / HLA-C eQTL. Columns: component, test_set, n_top, obs_overlap, matched_mean, matched_sd, fold_vs_matched, z_score, empirical_p (one-sided), empirical_p_reported (string). MAF bins [0, 0.01, 0.05, 0.10, 0.20, 0.30, 1.0]; distance bins [0, 25, 50, 100, 150, 200, 260] kb. *Legend note*: empirical *p* = 0 in the numeric column indicates that no permutation exceeded the observed statistic and is reported as **P < 1×10⁻⁴** (the 10,000-permutation floor) in the `empirical_p_reported` column.

**Table S3 — c5 axis eQTL direction across GTEx tissues.**
Source: `eqtl_axis_test.py`, `liver_replication.py`, `external_tissue_replication.py`. Per (analysis_block, component, tissue, target_gene): n_overlap, mean / median NES, n_pos / n_neg / sign_coherence (= |mean(sign(NES))|). Two analysis blocks (column `analysis_block`):

- *all_components_primary*: 8 NMF components × 5 target genes, fitted in **Whole_Blood** and **Liver**.
- *c5_cross_tissue_replication*: c5 only × 5 target genes, fitted in **Cells_EBV-transformed_lymphocytes**, **Lung**, and **Skin_Sun_Exposed_Lower_leg** to test reproducibility of the c5 direction in additional tissues.

  **Rows**: 87. The c5 rows in the cross-tissue block are *not* duplicates of the Whole_Blood / Liver primary block; they describe the same component evaluated in different tissues.

**Table S4 — LDA vs NMF c5-analogue enrichment.**
Source: `supp_1_lda_comparison.py`. Side-by-side fold enrichment for NMF c5 and LDA c5-analogue (Hungarian-matched) on 6 test sets (5 eQTL gene sets + C-flip).

**Table S5 — NMF rank cophenetic and reconstruction error.**
Source: `supp_2_cophenetic_k.py`. Per k ∈ {2, 3, ..., 12}: cophenetic_corr, dispersion, mean_recon_err, std_recon_err over 20 random-initialization seeds.

**Table S6 — NMF seed stability per component.**
Source: `supp_3_seed_stability.py`. Per component (c0..c7): mean / min / std of Hungarian-matched H-row Pearson correlation across 49 seeds vs seed 0 reference.

**Table S7 — Tag SNP → NMF component mapping (S7A long, S7B wide).**
Source: `tag_snp_component_map.py`, with biological-role annotations added by `make_supplementary_package.py`. Four manuscript-relevant tag SNPs (`rs2596542` anchor; `rs2395029` HLA-B*57:01 surrogate in HCP5; `rs2244546` Lange-2013 HCV-HCC tag in HCP5; `rs11509487` MICB upstream indel — the paper-1 anchor) are intersected with the 8 NMF components. The `in_NMF_set` column is True for `rs2395029` and `rs2244546`, which are biallelic SNVs that pass the MAF ≥ 0.05 filter and are direct rows of the H matrix. It is False for `rs2596542` (anchor, excluded from H by construction; used only to define the carrier set) and for `rs11509487` (a 9-bp insertion filtered from the biallelic SNV matrix). For these two, the H values reported are not for the tag itself but for the *nearest in-window NMF SNV* (`nearest_NMF_SNV_position` column); the column `distance_bp_to_nearest_NMF_SNV` makes this explicit (rs2596542: 374 bp; rs11509487: 43 bp). Two CSV variants are provided:

- `Table_S7A_tags_long.csv` — long format (1 row per (tag SNP × component); 32 rows). Columns: rsID, biological_role, grch37_pos, in_NMF_set, component, H_loading, rank_within_component_descending, rank_percentile_within_component, in_top_5pct_of_H_loadings, nearest_NMF_SNV_position, distance_bp_to_nearest_NMF_SNV, distance_bp_to_nearest_top5pct_SNV.
- `Table_S7B_tags_wide.csv` — wide format (1 row per tag SNP; 4 rows). Columns: rsID, biological_role, grch37_pos, in_NMF_set, H_c0..H_c7, in_top5pct_c0..c7, best_component, best_H_loading, n_components_in_top5pct.

Manuscript-relevant findings (best component = argmax H_loading; for tag SNPs not in the NMF matrix, H values are at the nearest in-NMF SNV — see column `distance_bp_to_nearest_NMF_SNV`):

- **rs2596542** (anchor; excluded from H by construction; nearest NMF SNV 374 bp away): nearest SNV's argmax is c3 (H = 0.33, rank 13 % within c3); H_c5 = 0.04 (rank 51 % within c5).
- **rs2395029** (HLA-B*57:01 tag; in_NMF_set = True): best component **c1** (H = 0.47, rank 32 % within c1); H_c5 = 0 (rank 64 % within c5; not in the top 5 % of c5).
- **rs2244546** (Lange 2013 tag; in_NMF_set = True): best component **c3** (H = 0.10, rank 59 % within c3), with smaller secondary loadings on c0 and c6; H_c5 = 0 (not in the top 5 % of c5).
- **rs11509487** (MICB 9-bp indel, paper-1 anchor; filtered as non-biallelic; nearest in-NMF SNV at 31,475,127, 43 bp away): nearest SNV's argmax is c2 (H = 0.49, **rank 12 % within c2**); H_c5 = 0.37 (rank 34 % within c5; not in the top 5 % of c5).

None of the four manuscript-relevant tag SNPs (or their nearest in-NMF SNVs) fall in the **top 5 %** of any component's H loading distribution, and in particular none of them fall in the top 5 % of **c5** (the HLA-B / HLA-C-variable axis). This supports the manuscript claim that the c5 axis is not captured by previously published HCV-HCC tag SNPs in the MICA / HCP5 / MICB region — c5 represents a haplotypic background that is uncharacterised in the existing tag-SNP literature.

**Table S8 — Cross-tissue NES direction analysis (Axis I + Axis II).**
Source: `axis_sign_coherence_matrix.py`. Companion to Figure 3 panel (c).

- `Table_S8_sign_coh.csv` — for every (component × tissue × gene) combination in {c4, c5, c6} × {Whole_Blood, Liver, Cells_EBV-transformed_lymphocytes (LCL), Lung, Skin_Sun_Exposed_Lower_leg, Colon_Transverse} × {MICA, HLA-B, HLA-C}, reports n_eqtl_in_window, n_top_overlap, n_pos, n_neg, sign_coherence, mean_NES, median_NES. Axis I cells (c4 / c6 × MICA): every one of 12 reported cells shows sign coherence = 1.00 in the negative direction (mean NES range −0.35 in Whole Blood to −0.72 in Colon Transverse; combined 2,731 SNV-tissue observations, n_pos = 0). Axis II cells (c5 × HLA-B, c5 × HLA-C): all 4–5 observed tissues per gene show sign coherence = 1.00 (HLA-B positive, HLA-C negative). Off-axis cells (e.g. c4 × HLA-C, c6 × HLA-C, c5 × MICA) included for completeness and show appreciably lower sign coherence (0.54–0.97), demonstrating that NMF cleanly separates the two axes by primary gene target rather than capturing a single hand-of-MHC haplotype.

---

## Excluded from supplementary (decisions log)

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
│   ├── figure_S4_comt_flipflop.{pdf,png}
│   └── _archive/                                     ← legacy figures (not cited from manuscript)
│       └── figure_2_branch_composition.{pdf,png}     ← old k = 2..5 hierarchical-clustering Figure 2
└── tables/
    ├── Table_S1_top_SNVs.csv
    ├── Table_S2_permutation.csv
    ├── Table_S3_GTEx_NES.csv
    ├── Table_S4_LDA_vs_NMF.csv
    ├── Table_S5_rank.csv
    ├── Table_S6_seed.csv
    ├── Table_S7A_tags_long.csv
    ├── Table_S7B_tags_wide.csv
    └── Table_S8_sign_coh.csv
```
