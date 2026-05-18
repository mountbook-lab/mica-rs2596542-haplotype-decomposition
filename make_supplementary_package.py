"""Assemble the AJHG-target supplementary package.

Outputs (under supplementary/):
  figures/figure_S1_lda_vs_nmf.{pdf,png}
  figures/figure_S2_cophenetic_k.{pdf,png}
  figures/figure_S3_seed_stability.{pdf,png}
  figures/figure_S4_liri_weighted.{pdf,png}
  figures/figure_S5_liri_hlaa_diagnostic.{pdf,png}

  tables/table_S1_regression_models.csv      ← assembled here (multi-section)
  tables/table_S2_component_top_snvs.csv     ← assembled here from nmf_H_26_k8.parquet
  tables/table_S3_matched_permutation.csv    ← copied from permutation_all_components.csv
  tables/table_S4_gtex_tissue_eqtl_NES.csv   ← merged Whole_Blood + Liver + LCL/Lung/Skin
  tables/table_S5_lda_vs_nmf.csv             ← copied from supp_1_lda_vs_nmf.csv
  tables/table_S6_nmf_rank_cophenetic.csv    ← copied from supp_2_cophenetic_k.csv
  tables/table_S7_nmf_seed_stability.csv     ← copied from supp_3_seed_stability_per_component.csv

  SUPPLEMENTARY_INDEX.md                      ← caption + source + key numbers per item

Survival (S6) is intentionally excluded — n = 16 events, dropped from manuscript on 2026-05-08.
"""

from __future__ import annotations

import os

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent
SRC = REPO / "results"
SUPP = REPO / "supplementary"
SUPP_FIG = SUPP / "figures"
SUPP_TAB = SUPP / "tables"

EXPR = str(Path(os.environ.get("LIRI_DATA_DIR", "data/liri_jp")) / "results/expression_matrix/tpm_gene_named.csv")
CLIN = str(Path(os.environ.get("LIRI_DATA_DIR", "data/liri_jp")) / "data/clinical_data/EGA_clinical_matched.csv")
LIRI_C5 = SRC / "liri_c5_score.csv"
SUPP4 = SRC / "supp_4_liri_weighted_per_sample.csv"
IMMUNE_GENES = ["PTPRC","CD3D","CD8A","IFNG","B2M","TAP1"]


def fit_table_s1():
    """Assemble Table S1 — regression models.

    Section A:  HLA-C ~ c5 unweighted, M0 → M3 progressive adjustment
    Section B:  HLA-C ~ c5 (full M3) — score variants {unweighted, weighted, NNLS}
    Section C:  Each gene (HLA-B / HLA-C / MICA / MICB / HCG27 / HCP5 / HLA-A) ~ c5 unweighted, full M3
    Section D:  VIF (c5_score_z, log2_HLA_A) for the full-model design matrix
    """
    df = pd.read_csv(LIRI_C5)
    clin = pd.read_csv(CLIN).rename(columns={"clinical_sample_id":"rk_id"})
    df = df.merge(clin[["rk_id","age","gender","viral_status","T_stage"]], on="rk_id")
    # Complete-case filter dropping 5 of 122 LIRI-JP donors:
    #   - 2 with `viral_status == "HBV, HCV"` (co-infection — not in {HBV/HCV/NBNC})
    #   - 3 with whitespace-suffix entries (`HBV `, `HCV `) that fall outside the strict isin filter
    # Rationale: the multivariable HLA-C model requires viral_status to be a
    # 3-level categorical for HCV / NBNC dummy coding; co-infection and
    # whitespace-irregular rows would either inflate categories or be silently
    # dropped by statsmodels. Documented in supplementary/SUPPLEMENTARY_INDEX.md
    # § "LIRI-JP analytic sample-size note (n = 122 vs n = 117)".
    df = df[df["T_stage"].notna() & df["age"].notna()
             & df["viral_status"].isin(["HBV","HCV","NBNC"])].copy()

    # immune markers
    expr = pd.read_csv(EXPR)
    gcol = expr.columns[0]; expr_samples = list(expr.columns[1:])
    feats = expr[gcol].astype(str)
    for g in IMMUNE_GENES + ["HLA-A","HLA-B","HLA-C","MICA","MICB","HCG27","HCP5"]:
        m = (feats == g)
        if not m.any(): continue
        row = expr.loc[m].iloc[0]
        df[g] = df["rk_id"].map(lambda rk: row[rk] if rk in expr_samples else np.nan)
        df[f"log2_{g}"] = np.log2(df[g].astype(float) + 1)

    cs = df["c5_score"].astype(float)
    df["c5_score_z"] = (cs - cs.mean()) / cs.std(ddof=0)
    df["sex_M"]       = (df["gender"] == "M").astype(int)
    df["viral_HCV"]   = (df["viral_status"] == "HCV").astype(int)
    df["viral_NBNC"]  = (df["viral_status"] == "NBNC").astype(int)
    df["age"]         = df["age"].astype(float)
    df["T_stage"]     = df["T_stage"].astype(float)

    # immune PC1
    Z = df[[f"log2_{g}" for g in IMMUNE_GENES]].copy()
    for c in Z.columns:
        Z[c] = (Z[c] - Z[c].mean()) / Z[c].std(ddof=0)
    Zv = Z.dropna().values
    U, S, Vt = np.linalg.svd(Zv, full_matrices=False)
    df["immune_PC1"] = (Z.values @ Vt[0])

    rows = []

    def fit(Y_col, X_cols, score_col, section, label):
        sub = df[[Y_col, score_col] + [c for c in X_cols if c != score_col]].dropna()
        # rename score col to "score_z" for unified tables (use as-is if already z)
        if score_col != "score_z":
            sub = sub.copy()
            s = sub[score_col].astype(float)
            sub["score_z"] = (s - s.mean()) / s.std(ddof=0)
        X = sub[["score_z"] + [c for c in X_cols if c != score_col]]
        m = sm.OLS(sub[Y_col], sm.add_constant(X)).fit()
        ci_lo, ci_hi = m.conf_int().loc["score_z"].tolist()
        return {
            "section": section,
            "label":   label,
            "outcome": Y_col,
            "score":   score_col,
            "n":       int(m.nobs),
            "beta":    float(m.params["score_z"]),
            "se":      float(m.bse["score_z"]),
            "t":       float(m.tvalues["score_z"]),
            "p":       float(m.pvalues["score_z"]),
            "ci_lo_95": float(ci_lo), "ci_hi_95": float(ci_hi),
            "r2":      float(m.rsquared),
            "covariates_in_model": ", ".join(X.columns.tolist()),
        }

    # ---------- Section A: HLA-C ~ c5 unweighted, M0..M3 ----------
    base = ["T_stage","viral_HCV","viral_NBNC","age","sex_M"]
    A = [
        fit("log2_HLA-C", ["c5_score_z"], "c5_score_z",
             "A", "M0: c5 only"),
        fit("log2_HLA-C", ["c5_score_z"] + base, "c5_score_z",
             "A", "M1: + clinical"),
        fit("log2_HLA-C", ["c5_score_z"] + base + ["immune_PC1"], "c5_score_z",
             "A", "M2: + immune_PC1"),
        fit("log2_HLA-C", ["c5_score_z"] + base + ["immune_PC1","log2_HLA-A"], "c5_score_z",
             "A", "M3: + immune_PC1 + HLA-A"),
    ]
    rows.extend(A)

    # ---------- Section B: HLA-C ~ c5 (M3 full) — score variants ----------
    # Pull weighted / NNLS from supp_4_liri_weighted_per_sample.csv
    if SUPP4.exists():
        s4 = pd.read_csv(SUPP4)
        # Bring the alternative scores in
        s4 = s4[["rk_id","c5_unweighted","c5_weighted","c5_nnls"]]
        df_b = df.merge(s4, on="rk_id", how="left")
        for sc, lab in [("c5_unweighted","unweighted (top-5%)"),
                          ("c5_weighted",  "H-loading weighted (all SNVs)"),
                          ("c5_nnls",      "NNLS-projected")]:
            sub = df_b[["log2_HLA-C", sc] + base + ["immune_PC1","log2_HLA-A"]].dropna()
            if len(sub) < 10: continue
            s = sub[sc].astype(float)
            sub = sub.copy()
            sub["score_z"] = (s - s.mean()) / s.std(ddof=0)
            X = sub[["score_z"] + base + ["immune_PC1","log2_HLA-A"]]
            m = sm.OLS(sub["log2_HLA-C"], sm.add_constant(X)).fit()
            ci_lo, ci_hi = m.conf_int().loc["score_z"].tolist()
            rows.append({
                "section": "B",
                "label":   f"M3 score variant: {lab}",
                "outcome": "log2_HLA-C",
                "score":   sc,
                "n":       int(m.nobs),
                "beta":    float(m.params["score_z"]),
                "se":      float(m.bse["score_z"]),
                "t":       float(m.tvalues["score_z"]),
                "p":       float(m.pvalues["score_z"]),
                "ci_lo_95": float(ci_lo), "ci_hi_95": float(ci_hi),
                "r2":      float(m.rsquared),
                "covariates_in_model": ", ".join(X.columns.tolist()),
            })

    # ---------- Section C: each gene as outcome, M3 full ----------
    full_covars = base + ["immune_PC1","log2_HLA-A"]
    for g in ["HLA-B","HLA-C","MICA","MICB","HCG27","HCP5","HLA-A"]:
        Y = f"log2_{g}"
        if Y not in df.columns: continue
        if g == "HLA-A":
            covars = base + ["immune_PC1"]   # avoid HLA-A on both sides
            label = "M3 (no HLA-A on RHS): outcome = HLA-A"
        else:
            covars = full_covars
            label = f"M3 full: outcome = {g}"
        rows.append(fit(Y, ["c5_score_z"] + covars, "c5_score_z", "C", label))

    # ---------- Section D: VIF for full M3 design ----------
    Xm3 = df[["c5_score_z"] + full_covars].dropna()
    def vif(Xdf, target):
        y = Xdf[target]; X = sm.add_constant(Xdf.drop(columns=[target]))
        m = sm.OLS(y, X).fit()
        return 1.0 / (1.0 - m.rsquared) if m.rsquared < 1 else float("inf")

    vif_rows = []
    for col in Xm3.columns:
        vif_rows.append({"section":"D", "label":"VIF (full M3 design)",
                          "predictor":col, "VIF": float(vif(Xm3, col))})
    df_vif = pd.DataFrame(vif_rows)

    df_main = pd.DataFrame(rows)
    df_main = df_main[["section","label","outcome","score","n",
                          "beta","se","t","p","ci_lo_95","ci_hi_95","r2",
                          "covariates_in_model"]]
    return df_main, df_vif


def extract_top_snvs():
    """Table S2 — top 5% SNVs per NMF component with annotations."""
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
    print("[Table S1] regression models ...")
    df_main, df_vif = fit_table_s1()
    df_main.to_csv(SUPP_TAB / "table_S1_regression_models.csv", index=False)
    df_vif.to_csv(SUPP_TAB / "table_S1d_VIF.csv", index=False)
    print(f"  wrote table_S1_regression_models.csv  ({len(df_main)} rows)")
    print(f"  wrote table_S1d_VIF.csv               ({len(df_vif)} rows)")

    print("[Table S2] component top SNVs ...")
    df_snv = extract_top_snvs()
    df_snv["H_loading"] = df_snv["H_loading"].round(4)
    df_snv.to_csv(SUPP_TAB / "table_S2_component_top_snvs.csv", index=False)
    print(f"  wrote table_S2_component_top_snvs.csv ({len(df_snv)} rows)")

    print("[Table S3] matched permutation ...")
    perm = pd.read_csv(SRC / "permutation_all_components.csv")
    # Replace empirical_p = 0 (10,000-perm floor) with the floor string.
    n_perm_floor = 1.0 / 10_000
    def _fmt_p(p):
        if pd.isna(p): return ""
        return f"<{n_perm_floor:.0e}" if p <= 0 else f"{p:.4f}"
    perm["empirical_p_reported"] = perm["empirical_p"].apply(_fmt_p)
    perm.to_csv(SUPP_TAB / "table_S3_matched_permutation.csv", index=False)
    n_floor = int((perm["empirical_p"] <= 0).sum())
    print(f"  wrote table_S3_matched_permutation.csv  ({len(perm)} rows; "
           f"{n_floor} rows hit the {n_perm_floor:.0e} floor)")

    print("[Table S4] GTEx tissue eQTL NES ...")
    df_tis = merge_gtex_tissue_NES()
    df_tis.to_csv(SUPP_TAB / "table_S4_gtex_tissue_eqtl_NES.csv", index=False)
    print(f"  wrote table_S4_gtex_tissue_eqtl_NES.csv ({len(df_tis)} rows)")

    print("[Table S5] LDA vs NMF ...")
    shutil.copy2(SRC / "supp_1_lda_vs_nmf.csv",
                  SUPP_TAB / "table_S5_lda_vs_nmf.csv")

    print("[Table S6] NMF rank cophenetic ...")
    shutil.copy2(SRC / "supp_2_cophenetic_k.csv",
                  SUPP_TAB / "table_S6_nmf_rank_cophenetic.csv")

    print("[Table S7] NMF seed stability ...")
    shutil.copy2(SRC / "supp_3_seed_stability_per_component.csv",
                  SUPP_TAB / "table_S7_nmf_seed_stability.csv")
    shutil.copy2(SRC / "supp_3_seed_stability.csv",
                  SUPP_TAB / "table_S7b_nmf_seed_stability_per_seed.csv")

    print("[Table S8] tag SNP → NMF component mapping ...")
    df_s8_long, df_s8_wide = build_table_s8()
    df_s8_long.to_csv(SUPP_TAB / "table_S8_tag_snp_component_mapping.csv",
                       index=False)
    df_s8_wide.to_csv(SUPP_TAB / "table_S8_tag_snp_component_mapping_wide.csv",
                       index=False)
    print(f"  wrote table_S8_tag_snp_component_mapping.csv "
           f"({len(df_s8_long)} long rows)")
    print(f"  wrote table_S8_tag_snp_component_mapping_wide.csv "
           f"({len(df_s8_wide)} tag SNPs)")

    print("[Table S9] LIRI sensitivity + multiple testing ...")
    if (SRC / "supp_6_liri_sensitivity.csv").exists():
        shutil.copy2(SRC / "supp_6_liri_sensitivity.csv",
                      SUPP_TAB / "table_S9_liri_sensitivity.csv")
        shutil.copy2(SRC / "supp_6_multiple_testing_section_C.csv",
                      SUPP_TAB / "table_S9b_multiple_testing.csv")
        print(f"  wrote table_S9_liri_sensitivity.csv")
        print(f"  wrote table_S9b_multiple_testing.csv")

    # ===== Figures =====
    print("\n[Figures] copy with S-prefix rename ...")
    fig_map = [
        ("figure_supp_1_lda_vs_nmf",        "figure_S1_lda_vs_nmf"),
        ("figure_supp_2_cophenetic_k",      "figure_S2_cophenetic_k"),
        ("figure_supp_3_seed_stability",    "figure_S3_seed_stability"),
        ("figure_supp_4_liri_weighted",     "figure_S4_liri_weighted"),
        ("figure_supp_5_liri_hlaa_conditioning",
                                              "figure_S5_liri_hlaa_diagnostic"),
        ("figure_supp_6_liri_sensitivity",  "figure_S6_liri_sensitivity"),
    ]
    for src_stem, dst_stem in fig_map:
        for ext in (".pdf", ".png"):
            shutil.copy2(SRC / f"{src_stem}{ext}",
                          SUPP_FIG / f"{dst_stem}{ext}")
        print(f"  {src_stem}.{{pdf,png}}  →  {dst_stem}.{{pdf,png}}")

    # ===== Index =====
    print("\n[INDEX] writing SUPPLEMENTARY_INDEX.md ...")
    write_index(df_main, df_snv, df_vif, df_tis)

    # ===== Summary =====
    n_fig = len(list(SUPP_FIG.glob("*.pdf")))
    n_tab = len(list(SUPP_TAB.glob("*.csv")))
    print(f"\nDone. supplementary/figures = {n_fig} PDFs, "
           f"supplementary/tables = {n_tab} CSVs")


def write_index(df_s1, df_s2, df_vif, df_s4):
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
    lines.append("## Supplementary Figures (7)")
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
    lines.append("**Figure S4 — LIRI-JP c5 score weighting sensitivity.**")
    lines.append("Source: `supp_4_liri_weighted.py`. "
                  "Three c5 score definitions are compared: unweighted top-5 % SNV burden "
                  "(n = 335 NMF SNVs available in LIRI), H-loading-weighted across all "
                  "7,116 NMF SNVs, and NNLS projection of LIRI dosage onto the 1000G "
                  "NMF basis. Panel (a) — unweighted vs weighted scatter; (b) unweighted "
                  "vs NNLS; (c) forest plot of c5 β on log₂(HLA-C TPM) under the full "
                  "M3 model (c5 + clinical + immune_PC1 + log₂(HLA-A); the same model "
                  "as Table S1 row B and Supplementary Figure S5). "
                  "**Key values** (n = 117; complete-case from 122 LIRI-JP donors with "
                  "matched genotype + RNA-seq, after dropping 5 donors with missing "
                  "T_stage or non-standard viral_status): Pearson r between unweighted "
                  "and weighted = 0.79; β (unweighted / weighted / NNLS) = "
                  "−0.193 / −0.160 / −0.175; all P < 0.005 in M3.")
    lines.append("")
    lines.append("**Figure S5 — HLA-A adjustment diagnostic for the LIRI-JP regression "
                  "of HLA-C on c5.**")
    lines.append("Source: `supp_5_liri_hlaa_conditioning.py`. "
                  "Diagnostic for whether HLA-A acts as a confounder of the c5 → HLA-C ↓ "
                  "association or as a precision-improving covariate capturing pan-class-I "
                  "shared expression variance. Panel (a) — c5 score vs log₂(HLA-A); "
                  "(b) — log₂(HLA-A) vs log₂(HLA-C); (c) — c5 score vs HLA-C residual "
                  "after regressing out HLA-A; (d) — forest plot of c5 β across M0–M3. "
                  "**Key values** (n = 117; complete-case from 122 LIRI-JP donors with "
                  "matched genotype + RNA-seq): c5 ⊥ HLA-A (Pearson r = −0.04, p = 0.64); "
                  "HLA-A explains 70 % of HLA-C variance; c5 vs HLA-C residual r = −0.37, "
                  "p = 5 × 10⁻⁵ (stronger than the raw r = −0.22); SE shrinkage M2 → M3 "
                  "= −30 %; β change = +19 % (suppression direction); "
                  "VIF(c5) = 1.10, VIF(HLA-A) = 2.10.")
    lines.append("")
    lines.append("**Figure S6 — LIRI-JP HLA-C regression sensitivity to analytic-sample "
                  "choice and multiple-testing context.**")
    lines.append("Source: `supp_6_liri_sensitivity.py`. "
                  "Three sensitivities: S1 (whitespace-trim of `viral_status` strings, "
                  "n: 117 → 120), S2 (S1 + re-assigning HBV+HCV co-infected donors to "
                  "HCV, n → 122), and S3 (c5 × viral_status interaction term on the "
                  "primary n = 117 set). Panel (a) — sensitivity forest plot; "
                  "(b) — −log₁₀(p) for the 7-outcome family (Table S1 section C) "
                  "under raw, BH-FDR, and Bonferroni adjustment. **Key values**: "
                  "primary β = −0.191 (p = 5×10⁻⁴), S1 β = −0.185 (p = 7×10⁻⁴), "
                  "S2 β = −0.187 (p = 5×10⁻⁴); c5 × HCV interaction term β = +0.18, "
                  "p = 0.24 (NS); c5 × NBNC term β = −0.05, p = 0.78 (NS). HLA-C is "
                  "the only outcome that survives both BH-FDR and Bonferroni "
                  "(adjusted p = 0.0037 by both methods, α = 0.05); the six other "
                  "tested gene outcomes are pure null after correction. "
                  "Locus specificity confirmed.")
    lines.append("")

    lines.append("**Figure S7 — Reference-panel replication and NMF extension of "
                  "the COMT Val158Met flip-flop example.**")
    lines.append("Source: `figure_S7_comt_flipflop.py` "
                  "(consumes `~/analyses/comt_flipflop_nmf/tables/comt_pair_metrics_by_pop.tsv`, "
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
    lines.append("## Supplementary Tables (10)")
    lines.append("")
    lines.append("**Table S1 — LIRI-JP regression models.**")
    lines.append("Built by `make_supplementary_package.py` from "
                  "`liri_c5_score.csv`, `EGA_clinical_matched.csv`, and "
                  "`tpm_gene_named.csv`. Multi-section table:")
    lines.append("")
    lines.append("- *Section A*: HLA-C ~ c5 (unweighted) progressive adjustment "
                  "M0 (c5 only) → M1 (+ clinical) → M2 (+ immune_PC1) → M3 (+ HLA-A).")
    lines.append("- *Section B*: HLA-C M3 full model with three c5 score variants "
                  "(unweighted / H-loading-weighted / NNLS-projected).")
    lines.append("- *Section C*: M3 full model fit to each MHC class-I / class-III "
                  "outcome (HLA-B, HLA-C, MICA, MICB, HCG27, HCP5, HLA-A).")
    f"  Reports columns: section, label, outcome, score, n, β, SE, t, p, "
    f"95 % CI, R², covariates_in_model."
    lines.append("- Companion `table_S1d_VIF.csv` reports VIF for every predictor in the "
                  "full M3 design matrix (c5 + clinical + immune_PC1 + HLA-A).")
    s1_n = len(df_s1)
    lines.append(f"  **Rows**: {s1_n} regression rows + 8 VIF rows.")
    lines.append("")
    lines.append("**Table S2 — Top-signature SNVs defining NMF components in rs2596542-T carrier haplotypes.**")
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
    lines.append("**Table S3 — Matched permutation for component-level enrichment.**")
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
    lines.append("**Table S4 — c5 axis eQTL direction across GTEx tissues.**")
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
    lines.append("**Table S5 — LDA vs NMF c5-analogue enrichment.**")
    lines.append("Source: `supp_1_lda_comparison.py`. "
                  "Side-by-side fold enrichment for NMF c5 and LDA c5-analogue (Hungarian-"
                  "matched) on 6 test sets (5 eQTL gene sets + C-flip).")
    lines.append("")
    lines.append("**Table S6 — NMF rank cophenetic and reconstruction error.**")
    lines.append("Source: `supp_2_cophenetic_k.py`. "
                  "Per k ∈ {2, 3, ..., 12}: cophenetic_corr, dispersion, "
                  "mean_recon_err, std_recon_err over 20 random-initialization seeds.")
    lines.append("")
    lines.append("**Table S7 — NMF seed stability per component.**")
    lines.append("Source: `supp_3_seed_stability.py`. "
                  "Per component (c0..c7): mean / min / std of Hungarian-matched H-row "
                  "Pearson correlation across 49 seeds vs seed 0 reference. "
                  "Companion `table_S7b_nmf_seed_stability_per_seed.csv` lists per-seed "
                  "matched correlations and MICA / HLA-B label recovery flags.")
    lines.append("")
    lines.append("**Table S8 — Tag SNP → NMF component mapping.**")
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
    lines.append("- `table_S8_tag_snp_component_mapping.csv` — long format "
                  "(1 row per (tag SNP × component); 32 rows). "
                  "Columns: rsID, biological_role, grch37_pos, in_NMF_set, "
                  "component, H_loading, rank_within_component_descending, "
                  "rank_percentile_within_component, in_top_5pct_of_H_loadings, "
                  "nearest_NMF_SNV_position, distance_bp_to_nearest_NMF_SNV, "
                  "distance_bp_to_nearest_top5pct_SNV.")
    lines.append("- `table_S8_tag_snp_component_mapping_wide.csv` — wide format "
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
    lines.append("**Table S9 — LIRI-JP regression sensitivities and multiple-testing.**")
    lines.append("Source: `supp_6_liri_sensitivity.py`. Two CSV files:")
    lines.append("")
    lines.append("- `table_S9_liri_sensitivity.csv` — per-row analysis label "
                  "(Primary n=117, S1 whitespace-trim n=120, S2 + HBV+HCV→HCV n=122, "
                  "S3 c5 × viral interaction terms), c5 β, SE, 95 % CI, p, R². "
                  "Confirms β stays in [−0.187, −0.191] across analytic-sample "
                  "definitions; interaction terms NS (HCV p = 0.24, NBNC p = 0.78).")
    lines.append("- `table_S9b_multiple_testing.csv` — Table S1 section C (7 gene "
                  "outcomes under M3 full model) augmented with BH-FDR and Bonferroni "
                  "adjusted p-values. HLA-C survives both corrections (adjusted p = "
                  "0.0037); the other six outcomes do not.")
    lines.append("")
    lines.append("**Table S10 — Cross-tissue NES direction analysis (Axis I + Axis II).**")
    lines.append("Source: `axis_sign_coherence_matrix.py`. Companion to Figure 3 panel (c).")
    lines.append("")
    lines.append("- `table_S10_axis_sign_coherence.csv` — for every "
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
    lines.append("## LIRI-JP analytic sample-size note (n = 122 vs n = 117)")
    lines.append("")
    lines.append("LIRI-JP included **122** donors with both germline genotype (rs2596542 "
                  "carriers + c5 top-signature SNVs available in `LIRI-JP_MHC_germline.vcf.gz`) "
                  "and tumour RNA-seq (EGA-derived TPM matrix; gene symbols indexed). "
                  "Univariate analyses (raw c5 score vs HLA-C, etc.) and the c5 tertile "
                  "boxplot in main Figure 5 (a) use n = 122.")
    lines.append("")
    lines.append("All **multivariable regression models** (Table S1 sections A–C, "
                  "Supplementary Figures S4 and S5) require complete data on c5 score, "
                  "expression, T_stage, viral_status (∈ {HBV, HCV, NBNC}), age, sex, "
                  "and HLA-A / immune-marker TPM. Five donors are dropped from the "
                  "complete-case subset because of missing T_stage or non-standard "
                  "viral_status entries (`HBV, HCV` co-infection, n = 2; whitespace-suffix "
                  "entries `HBV `, `HCV `, n = 3). The complete-case n is therefore **117**.")
    lines.append("")
    lines.append("Manuscript wording (insert into Methods):")
    lines.append("")
    lines.append("> *LIRI-JP included 122 donors with matched genotype and RNA-seq data; "
                  "complete-case multivariable regression models including all covariates "
                  "(T_stage, viral_status, age, sex, immune-marker expression, log₂(HLA-A)) "
                  "were performed in 117 donors after exclusion of 5 donors with missing "
                  "T_stage or non-standard viral_status entries.*")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Excluded from supplementary (decisions log)")
    lines.append("")
    lines.append("- **Survival analysis** (`liri_c5_survival.py`, "
                  "`results/liri_c5_survival_cox.csv`, "
                  "`results/figure_liri_c5_survival.{pdf,png}`) — DROPPED on 2026-05-08. "
                  "n = 117 with only 16 disease-specific death events, too underpowered "
                  "for primary AJHG submission. Files kept in repo for reproducibility "
                  "and for re-analysis if a replication cohort becomes available, but "
                  "not cited, summarised, or supplied as Supplementary in the manuscript.")
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
    lines.append("│   ├── figure_S4_liri_weighted.{pdf,png}")
    lines.append("│   ├── figure_S5_liri_hlaa_diagnostic.{pdf,png}")
    lines.append("│   ├── figure_S6_liri_sensitivity.{pdf,png}")
    lines.append("│   ├── figure_S7_comt_flipflop.{pdf,png}")
    lines.append("│   └── _archive/")
    lines.append("│       ├── figure_2_branch_composition.{pdf,png}")
    lines.append("│       └── figure_S_coarse_carrier_space_k2_k5.{pdf,png}")
    lines.append("└── tables/")
    lines.append("    ├── table_S1_regression_models.csv")
    lines.append("    ├── table_S1d_VIF.csv")
    lines.append("    ├── table_S2_component_top_snvs.csv")
    lines.append("    ├── table_S3_matched_permutation.csv")
    lines.append("    ├── table_S4_gtex_tissue_eqtl_NES.csv")
    lines.append("    ├── table_S5_lda_vs_nmf.csv")
    lines.append("    ├── table_S6_nmf_rank_cophenetic.csv")
    lines.append("    ├── table_S7_nmf_seed_stability.csv")
    lines.append("    ├── table_S7b_nmf_seed_stability_per_seed.csv")
    lines.append("    ├── table_S8_tag_snp_component_mapping.csv")
    lines.append("    ├── table_S8_tag_snp_component_mapping_wide.csv")
    lines.append("    ├── table_S9_liri_sensitivity.csv")
    lines.append("    ├── table_S9b_multiple_testing.csv")
    lines.append("    └── table_S10_axis_sign_coherence.csv")
    lines.append("```")
    lines.append("")

    text = "\n".join(lines)
    (SUPP / "SUPPLEMENTARY_INDEX.md").write_text(text)
    print(f"  wrote SUPPLEMENTARY_INDEX.md ({len(lines)} lines)")


if __name__ == "__main__":
    main()
