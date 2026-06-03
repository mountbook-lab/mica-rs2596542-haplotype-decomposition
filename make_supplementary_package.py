"""Assemble the AJHG-target supplementary package.

Outputs (under supplementary/):
  figures/figure_S1_lda_vs_nmf.{pdf,png}
  figures/figure_S2_cophenetic_k.{pdf,png}
  figures/figure_S3_seed_stability.{pdf,png}
  figures/figure_S4_comt_flipflop.{pdf,png}

  tables/  — final numbering (regression / VIF / per-seed / LIRI tables removed):
  Table_S1_top_SNVs.csv      ← assembled here from nmf_H_26_k8.parquet      (was table_S2)
  Table_S2_permutation.csv   ← permutation_all_components.csv               (was table_S3)
  Table_S3_GTEx_NES.csv      ← merged GTEx Whole_Blood/Liver/LCL/Lung/Skin  (was table_S4)
  Table_S4_LDA_vs_NMF.csv    ← supp_1_lda_vs_nmf.csv                        (was table_S5)
  Table_S5_rank.csv          ← supp_2_cophenetic_k.csv                      (was table_S6)
  Table_S6_seed.csv          ← supp_3_seed_stability_per_component.csv      (was table_S7)
  Table_S7A_tags_long.csv    ← tag SNP → component mapping, long            (was table_S8)
  Table_S7B_tags_wide.csv    ← tag SNP → component mapping, wide            (was table_S8 wide)
  Table_S8_sign_coh.csv      ← results/axis_sign_coherence_matrix.csv       (was table_S10)

  SUPPLEMENTARY_INDEX.md                      ← caption + source + key numbers per item
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
SRC = REPO / "results"
SUPP = REPO / "supplementary"
SUPP_FIG = SUPP / "figures"
SUPP_TAB = SUPP / "tables"




def extract_top_snvs():
    """Table S1 — top 5% SNVs per NMF component with annotations."""
    H = pd.read_parquet(SRC / "nmf_H_26_k8.parquet")
    nmf_pos = np.array([int(c) for c in H.columns])
    h_mat = H.values   # (8, 7116)
    K, n_snv = h_mat.shape
    top_n = max(1, int(0.05 * n_snv))

    # Annotations
    region = pd.read_csv(SRC / "region_per_partner_26.csv")
    region_idx = region.set_index("partner_pos")
    cflip = set(region.loc[region["any_C_flip"] == True,
                              "partner_pos"].astype(int))
    classflip = set(region.loc[region["any_class_flip"] == True,
                                  "partner_pos"].astype(int))

    eqtl = pd.read_csv(SRC / "eqtl_targets_per_gene.csv")
    eqtl_sets = {g: set(eqtl[eqtl["gene"] == g]["b37_pos"].astype(int))
                  for g in eqtl["gene"].unique()}

    # Approximate gene-region by position (GRCh37 chr6, well-known boundaries)
    def gene_region(p):
        if   31_366_595 <= p < 31_383_500: return "MICA"
        elif 31_462_660 <= p < 31_478_900: return "MICB"
        elif 31_321_649 <= p < 31_324_989: return "HLA-B"
        elif 31_236_526 <= p < 31_239_869: return "HLA-C"
        elif 31_165_537 <= p < 31_171_745: return "HCG27"
        elif 31_368_479 <= p < 31_445_283: return "HCP5"   # rough
        elif p < 31_165_537:               return "5'_of_HCG27"
        elif p > 31_478_900:               return "3'_of_MICB"
        else:                              return "intergenic_class_I"

    # Build long-format rows
    out = []
    for c in range(K):
        h = h_mat[c]
        order = np.argsort(h)[-top_n:][::-1]   # descending
        for rank, idx in enumerate(order, start=1):
            p = int(nmf_pos[idx])
            row = {
                "component": f"c{c}",
                "rank_within_component": rank,
                "position_grch37": p,
                "H_loading": float(h[idx]),
                "gene_region": gene_region(p),
                "C_flip":     bool(p in cflip),
                "class_flip": bool(p in classflip),
            }
            for g in ["MICA","MICB","HLA-B","HLA-C","HCG27","HCP5"]:
                row[f"eQTL_{g}"] = bool(p in eqtl_sets.get(g, set()))
            out.append(row)
    return pd.DataFrame(out)


def build_table_s8():
    """Tag SNP → NMF component mapping table.

    Returns (long_format_df, wide_format_df).

    Long format: 1 row per (tag SNP, component); columns include
    biological role, position, H_weight, rank, top-5 % indicator,
    nearest-top-5 % distance.

    Wide format: 1 row per tag SNP; columns include H_c0..H_c7 weights,
    in_top5pct_c0..c7 flags, best_component, n_components_top5pct.
    """
    src = pd.read_csv(SRC / "tag_snp_component_map.csv")

    # Map each tag-SNP label to a clean rsID + biological role
    # (the existing CSV labels embed the role parenthetically)
    role_map = {
        "rs2596542 (anchor)": (
            "rs2596542",
            "MICA upstream (5′ flanking; ~2 kb of MICA); "
            "Kumar 2011 HCV-HCC GWAS hit; *anchor of NMF carrier set, "
            "excluded from H matrix*"),
        "rs2395029 (HCP5, B*57:01 tag)": (
            "rs2395029",
            "HCP5 intronic; HLA-B*57:01 surrogate; PSC, drug-hypersensitivity, "
            "HIV-control reported"),
        "rs2244546 (Lange 2013 tag)": (
            "rs2244546",
            "HCP5 region; Lange 2013 HCV-HCC tag with consistent allele direction "
            "across European and Japanese cohorts"),
        "rs11509487 (MICB, paper 1)": (
            "rs11509487",
            "MICB upstream (paper-1 anchor); 9-bp insertion (chr6:31,475,084 "
            "C>CTGGGGTGA); class-I-like haplotype tag"),
    }
    src["rsID"]            = src["tag_snp"].map(lambda x: role_map[x][0])
    src["biological_role"] = src["tag_snp"].map(lambda x: role_map[x][1])

    long_df = src.rename(columns={
        "tag_pos": "grch37_pos",
        "H_weight": "H_loading",
        "rank_in_component": "rank_within_component_descending",
        "rank_percentile":   "rank_percentile_within_component",
        "in_top5pct":        "in_top_5pct_of_H_loadings",
        "nearest_top5pct_dist_bp": "distance_bp_to_nearest_top5pct_SNV",
        "nearest_NMF_pos":         "nearest_NMF_SNV_position",
        "nearest_NMF_dist_bp":     "distance_bp_to_nearest_NMF_SNV",
    })
    long_df = long_df[[
        "rsID", "biological_role", "grch37_pos", "in_NMF_set",
        "component",
        "H_loading", "rank_within_component_descending",
        "rank_percentile_within_component", "in_top_5pct_of_H_loadings",
        "nearest_NMF_SNV_position",
        "distance_bp_to_nearest_NMF_SNV",
        "distance_bp_to_nearest_top5pct_SNV",
    ]]

    # Wide pivot: 1 row per tag SNP
    wide_rows = []
    for rs in long_df["rsID"].drop_duplicates():
        sub = long_df[long_df["rsID"] == rs]
        first = sub.iloc[0]
        row = {
            "rsID":              rs,
            "biological_role":   first["biological_role"],
            "grch37_pos":        int(first["grch37_pos"]),
            "in_NMF_set":        bool(first["in_NMF_set"]),
        }
        # Per-component H_loading + top-5 % flag
        sub_idx = sub.set_index("component")
        for c in range(8):
            if c in sub_idx.index:
                h = float(sub_idx.loc[c, "H_loading"])
                in5 = bool(sub_idx.loc[c, "in_top_5pct_of_H_loadings"])
            else:
                h, in5 = float("nan"), False
            row[f"H_c{c}"] = h
            row[f"in_top5pct_c{c}"] = in5

        # Summary
        H_vec = np.array([row[f"H_c{c}"] for c in range(8)])
        if np.all(np.isnan(H_vec)):
            row["best_component"]            = ""
            row["best_H_loading"]             = float("nan")
            row["n_components_in_top5pct"]    = 0
        else:
            best = int(np.nanargmax(H_vec))
            row["best_component"]            = f"c{best}"
            row["best_H_loading"]             = float(H_vec[best])
            row["n_components_in_top5pct"]    = int(sum(
                row[f"in_top5pct_c{c}"] for c in range(8)))
        wide_rows.append(row)
    wide_df = pd.DataFrame(wide_rows)
    return long_df, wide_df


def merge_gtex_tissue_NES():
    """Table S4 — c5 axis NES across GTEx tissues.

    Two analysis blocks (controlled by the `analysis_block` column to
    avoid the apparent c5 row duplication between Whole_Blood / Liver
    primary and the cross-tissue replication block):

      `all_components_primary`         : 8 NMF components × 5 target genes,
                                          fitted in Whole_Blood and Liver only.
      `c5_cross_tissue_replication`    : c5 only × 5 target genes, fitted in
                                          Cells_EBV-transformed_lymphocytes,
                                          Lung, and Skin_Sun_Exposed_Lower_leg.
    Rows in the cross-tissue block are *not* duplicates of Whole_Blood / Liver
    rows; they describe the same c5 component evaluated in different
    tissues to test reproducibility of the c5 → HLA-B / HLA-C / HCG27 axis.
    """
    blood = pd.read_csv(SRC / "eqtl_NES_per_component.csv")
    blood = blood.rename(columns={"target_gene":"gene"})
    blood["tissue"] = "Whole_Blood"
    blood["analysis_block"] = "all_components_primary"
    blood = blood[["analysis_block","component","tissue","gene","n_overlap",
                     "mean_NES","median_NES","n_pos_NES","n_neg_NES",
                     "sign_coherence"]
                    ].rename(columns={"n_pos_NES":"n_pos","n_neg_NES":"n_neg"})

    liver = pd.read_csv(SRC / "liver_NES_per_component.csv")
    if "n_eqtl_in_window" in liver.columns:
        liver = liver.drop(columns=["n_eqtl_in_window"])
    liver["analysis_block"] = "all_components_primary"
    liver = liver[["analysis_block","component","tissue","gene","n_overlap",
                     "mean_NES","median_NES","n_pos","n_neg","sign_coherence"]]

    other = pd.read_csv(SRC / "c5_external_tissue_NES.csv")
    other = other.rename(columns={"n_overlap_c5":"n_overlap"})
    # Drop Whole_Blood rows from the cross-tissue file to avoid duplicating
    # the primary-block c5 rows already represented above.
    other = other[other["tissue"] != "Whole_Blood"].copy()
    other["component"] = 5
    other["analysis_block"] = "c5_cross_tissue_replication"
    other = other[["analysis_block","component","tissue","gene","n_overlap",
                     "mean_NES","median_NES","n_pos","n_neg","sign_coherence"]]

    merged = pd.concat([blood, liver, other], ignore_index=True, sort=False)
    return merged


def main():
    sys.stdout.reconfigure(line_buffering=True)
    SUPP_FIG.mkdir(parents=True, exist_ok=True)
    SUPP_TAB.mkdir(parents=True, exist_ok=True)

    # ===== Tables =====
    print("[Table S1] component top SNVs ...")
    df_snv = extract_top_snvs()
    df_snv["H_loading"] = df_snv["H_loading"].round(4)
    df_snv.to_csv(SUPP_TAB / "Table_S1_top_SNVs.csv", index=False)
    print(f"  wrote Table_S1_top_SNVs.csv ({len(df_snv)} rows)")

    print("[Table S2] matched permutation ...")
    perm = pd.read_csv(SRC / "permutation_all_components.csv")
    # Replace empirical_p = 0 (10,000-perm floor) with the floor string.
    n_perm_floor = 1.0 / 10_000
    def _fmt_p(p):
        if pd.isna(p): return ""
        return f"<{n_perm_floor:.0e}" if p <= 0 else f"{p:.4f}"
    perm["empirical_p_reported"] = perm["empirical_p"].apply(_fmt_p)
    perm.to_csv(SUPP_TAB / "Table_S2_permutation.csv", index=False)
    n_floor = int((perm["empirical_p"] <= 0).sum())
    print(f"  wrote Table_S2_permutation.csv  ({len(perm)} rows; "
           f"{n_floor} rows hit the {n_perm_floor:.0e} floor)")

    print("[Table S3] GTEx tissue eQTL NES ...")
    df_tis = merge_gtex_tissue_NES()
    df_tis.to_csv(SUPP_TAB / "Table_S3_GTEx_NES.csv", index=False)
    print(f"  wrote Table_S3_GTEx_NES.csv ({len(df_tis)} rows)")

    print("[Table S4] LDA vs NMF ...")
    shutil.copy2(SRC / "supp_1_lda_vs_nmf.csv",
                  SUPP_TAB / "Table_S4_LDA_vs_NMF.csv")

    print("[Table S5] NMF rank cophenetic ...")
    shutil.copy2(SRC / "supp_2_cophenetic_k.csv",
                  SUPP_TAB / "Table_S5_rank.csv")

    print("[Table S6] NMF seed stability ...")
    shutil.copy2(SRC / "supp_3_seed_stability_per_component.csv",
                  SUPP_TAB / "Table_S6_seed.csv")

    print("[Table S7A/S7B] tag SNP → NMF component mapping ...")
    df_s8_long, df_s8_wide = build_table_s8()
    df_s8_long.to_csv(SUPP_TAB / "Table_S7A_tags_long.csv", index=False)
    df_s8_wide.to_csv(SUPP_TAB / "Table_S7B_tags_wide.csv", index=False)
    print(f"  wrote Table_S7A_tags_long.csv ({len(df_s8_long)} long rows)")
    print(f"  wrote Table_S7B_tags_wide.csv ({len(df_s8_wide)} tag SNPs)")

    print("[Table S8] cross-tissue NES sign coherence ...")
    shutil.copy2(SRC / "axis_sign_coherence_matrix.csv",
                  SUPP_TAB / "Table_S8_sign_coh.csv")
    print("  wrote Table_S8_sign_coh.csv")

    # ===== Figures =====
    print("\n[Figures] copy with S-prefix rename ...")
    fig_map = [
        ("figure_supp_1_lda_vs_nmf",        "figure_S1_lda_vs_nmf"),
        ("figure_supp_2_cophenetic_k",      "figure_S2_cophenetic_k"),
        ("figure_supp_3_seed_stability",    "figure_S3_seed_stability"),
    ]
    for src_stem, dst_stem in fig_map:
        for ext in (".pdf", ".png"):
            shutil.copy2(SRC / f"{src_stem}{ext}",
                          SUPP_FIG / f"{dst_stem}{ext}")
        print(f"  {src_stem}.{{pdf,png}}  →  {dst_stem}.{{pdf,png}}")

    # ===== Index =====
    print("\n[INDEX] writing SUPPLEMENTARY_INDEX.md ...")
    write_index(df_snv, df_tis)

    # ===== Summary =====
    n_fig = len(list(SUPP_FIG.glob("*.pdf")))
    n_tab = len(list(SUPP_TAB.glob("*.csv")))
    print(f"\nDone. supplementary/figures = {n_fig} PDFs, "
           f"supplementary/tables = {n_tab} CSVs")


def write_index(df_s2, df_s4):
    """Write the SUPPLEMENTARY_INDEX.md with caption + key numbers per item."""
    lines = []
    lines.append("# Supplementary materials index — pre-AJHG submission")
    lines.append("")
    lines.append("Date: 2026-05-08.  Project: `mica-rs2596542-haplotype-decomposition`.")
    lines.append("")
    lines.append("All numbers in the captions below are taken directly from the corresponding "
                  "table CSV. Source scripts are in the parent repo.")
    lines.append("")
    lines.append("---")
    lines.append("")

    # Figures
    lines.append("## Supplementary Figures (4)")
    lines.append("")
    lines.append("**Figure S1 — Comparison of NMF and LDA decomposition on the same "
                  "rs2596542-T carrier haplotype matrix.**")
    lines.append("Source: `supp_1_lda_comparison.py`. "
                  "k = 8 LatentDirichletAllocation fitted on the same 2,111 × 7,116 "
                  "haplotype matrix; Hungarian matching of LDA components to NMF components "
                  "via H-row Pearson correlation. Panel (a) — matched correlation per "
                  "NMF component; (b) — LDA-c5 analogue eQTL / C-flip enrichment "
                  "alongside NMF-c5 (c5-like enrichment pattern is *partially retained* "
                  "under LDA). **Key values**: matched corr c4 / c5 / c6 = "
                  "0.99 / 0.87 / 0.96; LDA-c5 preserves MICA depletion (fold = 0) "
                  "and C-flip enrichment (fold = 1.39) but attenuates the localized "
                  "HLA-B *cis*-eQTL signal (NMF fold 2.14 → LDA fold 0.66). LDA "
                  "recovered the major component structure but smoothed the localized "
                  "HLA-B regulatory signature, supporting NMF as the appropriate method "
                  "for sparse, localized component-defining SNV signatures.")
    lines.append("")
    lines.append("**Figure S2 — NMF rank selection by consensus cophenetic correlation.**")
    lines.append("Source: `supp_2_cophenetic_k.py`. "
                  "k = 2..12 with 20 random initializations per rank; consensus matrix "
                  "C[i,j] = (fraction of seeds where haplotypes i and j share dominant "
                  "component); cophenetic correlation = Pearson r between cophenetic "
                  "distances on the average-linkage dendrogram and (1 − C). "
                  "Panel (a) — cophenetic correlation and Brunet-2004 dispersion vs k; "
                  "(b) — mean reconstruction error ‖X − WH‖_F. "
                  "**Key values**: cophenetic max = 0.998 at k = 2; k = 8 cophenetic = "
                  "0.963 (within 0.05 of max).")
    lines.append("")
    lines.append("**Figure S3 — NMF seed reproducibility at k = 8.**")
    lines.append("Source: `supp_3_seed_stability.py`. "
                  "50 independent NMF runs at k = 8 (init = random). Hungarian-matched H-row "
                  "Pearson correlations (anchored to seed 0); recovery of MICA-enriched "
                  "(fold > 3) and HLA-B-enriched (fold > 1.5) component labels per seed. "
                  "Panel (a) — boxplot of matched corr per component; (b) — MICA / HLA-B "
                  "fold across 49 matched seeds. "
                  "**Key values**: c4 / c5 / c6 mean matched corr = 0.98 / 0.96 / 0.99; "
                  "MICA-enriched recovered in 46/49 (94 %) seeds; "
                  "HLA-B-enriched recovered in 49/49 (100 %) seeds.")
    lines.append("")
    lines.append("**Figure S4 — Reference-panel replication and NMF extension of "
                  "the COMT Val158Met flip-flop example.**")
    lines.append("Source: `figure_S4_comt_flipflop.py` "
                  "(consumes `data/comt_flipflop_nmf/tables/comt_pair_metrics_by_pop.tsv`, "
                  "`comt_component_mixture_by_rs4680_status.tsv`, "
                  "`comt_component_mixture_by_pop.tsv`). "
                  "(A) Population-conditional signed LD and asymmetry component C for "
                  "rs4680 (Val158Met) and rs2097603 across 1000 Genomes populations. "
                  "Blue and brown indicate positive and negative signed values, respectively. "
                  "(B) Row-normalized k = 4 NMF component mixture conditional on rs4680 "
                  "allele status, with V indicating Val/REF carriers and M indicating "
                  "Met/ALT carriers. "
                  "(C) Population-level COMT NMF component mixture across all 5,008 haplotypes. "
                  "Compared with the MHC rs2596542 locus, the COMT decomposition is simpler "
                  "and more population-like, supporting broader applicability of carrier-set "
                  "decomposition while illustrating locus-dependent differences in functional "
                  "interpretability. COMT is included as an external reference-panel "
                  "proof-of-concept and is not cited in the main Results.")
    lines.append("")

    # Tables
    lines.append("## Supplementary Tables (9)")
    lines.append("")
    lines.append("**Table S1 — Top-signature SNVs defining NMF components in rs2596542-T carrier haplotypes.**")
    lines.append("Built by `make_supplementary_package.py` from "
                  "`results/nmf_H_26_k8.parquet` and "
                  "`results/region_per_partner_26.csv`. "
                  "Top-signature SNVs were defined as variants with H-matrix loadings in "
                  "the top 5 % of each component's H[c, ·] distribution (i.e., at or above "
                  "the 95th percentile), yielding 355 SNVs per component. The anchor SNV "
                  "rs2596542 was removed from the NMF input matrix after carrier selection "
                  "because all retained haplotypes carry the rs2596542-T allele at this site, "
                  "leaving it with zero variance in the analytic matrix. Columns: "
                  "`component` (c0..c7); `rank_within_component` (1 = highest H_loading, "
                  "descending); `position_grch37` (chr6 1000G phase 3 coordinate); "
                  "`H_loading` (rounded to 4 decimal places); `gene_region` "
                  "(MICA / MICB / HLA-B / HLA-C / HCG27 / HCP5 / 3'_of_MICB / "
                  "5'_of_HCG27 / intergenic_class_I); `C_flip`, `class_flip` "
                  "(boolean LD-reversal flags); `eQTL_<gene>` (boolean indicating "
                  "whether the SNV is a GTEx Whole_Blood *cis*-eQTL for that gene). "
                  f"**Rows**: {len(df_s2)} (8 components × 355).")
    lines.append("")
    lines.append("**Table S2 — Matched permutation for component-level enrichment.**")
    lines.append("Source: `permutation_all_components.py`. "
                  "Per (NMF component × test set) MAF×distance matched permutation, "
                  "10,000 perms. Test sets: r-flip / C-flip / class-flip / "
                  "MICA / MICB / HLA-B / HLA-C eQTL. "
                  "Columns: component, test_set, n_top, obs_overlap, matched_mean, "
                  "matched_sd, fold_vs_matched, z_score, empirical_p (one-sided), "
                  "empirical_p_reported (string). "
                  "MAF bins [0, 0.01, 0.05, 0.10, 0.20, 0.30, 1.0]; "
                  "distance bins [0, 25, 50, 100, 150, 200, 260] kb. "
                  "*Legend note*: empirical *p* = 0 in the numeric column indicates "
                  "that no permutation exceeded the observed statistic and is reported "
                  "as **P < 1×10⁻⁴** (the 10,000-permutation floor) in the "
                  "`empirical_p_reported` column.")
    lines.append("")
    lines.append("**Table S3 — c5 axis eQTL direction across GTEx tissues.**")
    lines.append("Source: `eqtl_axis_test.py`, `liver_replication.py`, "
                  "`external_tissue_replication.py`. "
                  "Per (analysis_block, component, tissue, target_gene): n_overlap, "
                  "mean / median NES, n_pos / n_neg / sign_coherence "
                  "(= |mean(sign(NES))|). Two analysis blocks (column `analysis_block`):")
    lines.append("")
    lines.append("- *all_components_primary*: 8 NMF components × 5 target genes, "
                  "fitted in **Whole_Blood** and **Liver**.")
    lines.append("- *c5_cross_tissue_replication*: c5 only × 5 target genes, fitted in "
                  "**Cells_EBV-transformed_lymphocytes**, **Lung**, and "
                  "**Skin_Sun_Exposed_Lower_leg** to test reproducibility of the c5 "
                  "direction in additional tissues.")
    lines.append("")
    lines.append(f"  **Rows**: {len(df_s4)}. The c5 rows in the cross-tissue block are "
                  "*not* duplicates of the Whole_Blood / Liver primary block; they "
                  "describe the same component evaluated in different tissues.")
    lines.append("")
    lines.append("**Table S4 — LDA vs NMF c5-analogue enrichment.**")
    lines.append("Source: `supp_1_lda_comparison.py`. "
                  "Side-by-side fold enrichment for NMF c5 and LDA c5-analogue (Hungarian-"
                  "matched) on 6 test sets (5 eQTL gene sets + C-flip).")
    lines.append("")
    lines.append("**Table S5 — NMF rank cophenetic and reconstruction error.**")
    lines.append("Source: `supp_2_cophenetic_k.py`. "
                  "Per k ∈ {2, 3, ..., 12}: cophenetic_corr, dispersion, "
                  "mean_recon_err, std_recon_err over 20 random-initialization seeds.")
    lines.append("")
    lines.append("**Table S6 — NMF seed stability per component.**")
    lines.append("Source: `supp_3_seed_stability.py`. "
                  "Per component (c0..c7): mean / min / std of Hungarian-matched H-row "
                  "Pearson correlation across 49 seeds vs seed 0 reference.")
    lines.append("")
    lines.append("**Table S7 — Tag SNP → NMF component mapping (S7A long, S7B wide).**")
    lines.append("Source: `tag_snp_component_map.py`, with biological-role annotations "
                  "added by `make_supplementary_package.py`. "
                  "Four manuscript-relevant tag SNPs (`rs2596542` anchor; "
                  "`rs2395029` HLA-B*57:01 surrogate in HCP5; `rs2244546` Lange-2013 "
                  "HCV-HCC tag in HCP5; `rs11509487` MICB upstream indel — the paper-1 "
                  "anchor) are intersected with the 8 NMF components. The `in_NMF_set` "
                  "column is True for `rs2395029` and `rs2244546`, which are biallelic "
                  "SNVs that pass the MAF ≥ 0.05 filter and are direct rows of the "
                  "H matrix. It is False for `rs2596542` (anchor, excluded from H by "
                  "construction; used only to define the carrier set) and for "
                  "`rs11509487` (a 9-bp insertion filtered from the biallelic SNV "
                  "matrix). For these two, the H values reported are not for the tag "
                  "itself but for the *nearest in-window NMF SNV* "
                  "(`nearest_NMF_SNV_position` column); the column "
                  "`distance_bp_to_nearest_NMF_SNV` makes this explicit "
                  "(rs2596542: 374 bp; rs11509487: 43 bp). "
                  "Two CSV variants are provided:")
    lines.append("")
    lines.append("- `Table_S7A_tags_long.csv` — long format "
                  "(1 row per (tag SNP × component); 32 rows). "
                  "Columns: rsID, biological_role, grch37_pos, in_NMF_set, "
                  "component, H_loading, rank_within_component_descending, "
                  "rank_percentile_within_component, in_top_5pct_of_H_loadings, "
                  "nearest_NMF_SNV_position, distance_bp_to_nearest_NMF_SNV, "
                  "distance_bp_to_nearest_top5pct_SNV.")
    lines.append("- `Table_S7B_tags_wide.csv` — wide format "
                  "(1 row per tag SNP; 4 rows). "
                  "Columns: rsID, biological_role, grch37_pos, in_NMF_set, "
                  "H_c0..H_c7, in_top5pct_c0..c7, best_component, "
                  "best_H_loading, n_components_in_top5pct.")
    lines.append("")
    lines.append("Manuscript-relevant findings (best component = argmax H_loading; "
                  "for tag SNPs not in the NMF matrix, H values are at the nearest "
                  "in-NMF SNV — see column `distance_bp_to_nearest_NMF_SNV`):")
    lines.append("")
    lines.append("- **rs2596542** (anchor; excluded from H by construction; nearest "
                  "NMF SNV 374 bp away): nearest SNV's argmax is c3 (H = 0.33, rank "
                  "13 % within c3); H_c5 = 0.04 (rank 51 % within c5).")
    lines.append("- **rs2395029** (HLA-B*57:01 tag; in_NMF_set = True): best "
                  "component **c1** (H = 0.47, rank 32 % within c1); H_c5 = 0 (rank "
                  "64 % within c5; not in the top 5 % of c5).")
    lines.append("- **rs2244546** (Lange 2013 tag; in_NMF_set = True): best component "
                  "**c3** (H = 0.10, rank 59 % within c3), with smaller secondary "
                  "loadings on c0 and c6; H_c5 = 0 (not in the top 5 % of c5).")
    lines.append("- **rs11509487** (MICB 9-bp indel, paper-1 anchor; filtered as "
                  "non-biallelic; nearest in-NMF SNV at 31,475,127, 43 bp away): "
                  "nearest SNV's argmax is c2 (H = 0.49, **rank 12 % within c2**); "
                  "H_c5 = 0.37 (rank 34 % within c5; not in the top 5 % of c5).")
    lines.append("")
    lines.append("None of the four manuscript-relevant tag SNPs (or their nearest "
                  "in-NMF SNVs) fall in the **top 5 %** of any component's H loading "
                  "distribution, and in particular none of them fall in the top 5 % "
                  "of **c5** (the HLA-B / HLA-C-variable axis). This supports the "
                  "manuscript claim that the c5 axis is not captured by previously "
                  "published HCV-HCC tag SNPs in the MICA / HCP5 / MICB region — c5 "
                  "represents a haplotypic background that is uncharacterised in the "
                  "existing tag-SNP literature.")
    lines.append("")
    lines.append("**Table S8 — Cross-tissue NES direction analysis (Axis I + Axis II).**")
    lines.append("Source: `axis_sign_coherence_matrix.py`. Companion to Figure 3 panel (c).")
    lines.append("")
    lines.append("- `Table_S8_sign_coh.csv` — for every "
                  "(component × tissue × gene) combination in "
                  "{c4, c5, c6} × {Whole_Blood, Liver, Cells_EBV-transformed_lymphocytes (LCL), "
                  "Lung, Skin_Sun_Exposed_Lower_leg, Colon_Transverse} × {MICA, HLA-B, HLA-C}, "
                  "reports n_eqtl_in_window, n_top_overlap, n_pos, n_neg, sign_coherence, "
                  "mean_NES, median_NES. Axis I cells (c4 / c6 × MICA): every one of "
                  "12 reported cells shows sign coherence = 1.00 in the negative direction "
                  "(mean NES range −0.35 in Whole Blood to −0.72 in Colon Transverse; "
                  "combined 2,731 SNV-tissue observations, n_pos = 0). Axis II cells "
                  "(c5 × HLA-B, c5 × HLA-C): all 4–5 observed tissues per gene show "
                  "sign coherence = 1.00 (HLA-B positive, HLA-C negative). Off-axis cells "
                  "(e.g. c4 × HLA-C, c6 × HLA-C, c5 × MICA) included for completeness and "
                  "show appreciably lower sign coherence (0.54–0.97), demonstrating that "
                  "NMF cleanly separates the two axes by primary gene target rather than "
                  "capturing a single hand-of-MHC haplotype.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Excluded from supplementary (decisions log)")
    lines.append("")
    lines.append("- **Internal PASS/FAIL summary** (`results/supp_combined_summary.csv`) — "
                  "internal QA artifact, not for publication. The thresholds it encodes "
                  "(matched corr > 0.7, cophenetic-within-0.05, score corr > 0.8) are "
                  "post-hoc rules-of-thumb chosen during analysis design and should not "
                  "be presented as pre-registered tests.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## File-system layout")
    lines.append("")
    lines.append("```")
    lines.append("supplementary/")
    lines.append("├── SUPPLEMENTARY_INDEX.md")
    lines.append("├── figures/")
    lines.append("│   ├── figure_S1_lda_vs_nmf.{pdf,png}")
    lines.append("│   ├── figure_S2_cophenetic_k.{pdf,png}")
    lines.append("│   ├── figure_S3_seed_stability.{pdf,png}")
    lines.append("│   ├── figure_S4_comt_flipflop.{pdf,png}")
    lines.append("│   └── _archive/")
    lines.append("│       ├── figure_2_branch_composition.{pdf,png}")
    lines.append("│       └── figure_S_coarse_carrier_space_k2_k5.{pdf,png}")
    lines.append("└── tables/")
    lines.append("    ├── Table_S1_top_SNVs.csv")
    lines.append("    ├── Table_S2_permutation.csv")
    lines.append("    ├── Table_S3_GTEx_NES.csv")
    lines.append("    ├── Table_S4_LDA_vs_NMF.csv")
    lines.append("    ├── Table_S5_rank.csv")
    lines.append("    ├── Table_S6_seed.csv")
    lines.append("    ├── Table_S7A_tags_long.csv")
    lines.append("    ├── Table_S7B_tags_wide.csv")
    lines.append("    └── Table_S8_sign_coh.csv")
    lines.append("```")
    lines.append("")

    text = "\n".join(lines)
    (SUPP / "SUPPLEMENTARY_INDEX.md").write_text(text)
    print(f"  wrote SUPPLEMENTARY_INDEX.md ({len(lines)} lines)")


if __name__ == "__main__":
    main()
