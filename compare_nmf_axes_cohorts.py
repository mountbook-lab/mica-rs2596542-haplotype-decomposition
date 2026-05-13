"""Compare NMF-axis proportions between 1000G EAS and LIRI-JP HCC.

NMF axes are defined on rs2596542-T carrier haplotypes (k=8, 26-pop;
results/nmf_H_26_k8.parquet).  This script projects haplotypes /
samples from independent cohorts onto the same H signature and
asks whether the *proportion* of each axis differs between

  - 1000 Genomes EAS 5 sub-pops (CHB / JPT / CHS / CDX / KHV)
  - LIRI-JP HCC patient cohort

separately for

  - Axis I  (MICA-stable / MICA↓ axis)  = c4 + c6
  - Axis II (HLA-B / HLA-C-variable; signed-LD reversal axis) = c5

Proportions are computed both
  (a) per haplotype  : 1000G uses phased haps directly; LIRI uses
      diploid dosage split halfway between the two haps before
      NNLS projection.
  (b) per diploid sample : diploid alt-dosage projected onto H by NNLS.

For each haplotype / sample we also record rs2596542 carrier status
so we can compare T-carrier vs non-carrier within each cohort, and
LIRI-side is further stratified by viral_status (HBV / HCV / NBNC).

Outputs:
  results/cohort_axis_proportions_perhap.csv
  results/cohort_axis_proportions_persample.csv
  results/cohort_axis_proportions_stats.csv
  results/figure_cohort_axis_proportions.{pdf,png}
"""

from __future__ import annotations

import gzip
import os
import subprocess
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import mannwhitneyu

REPO = Path(__file__).resolve().parent
H_PARQUET = REPO / "results/nmf_H_26_k8.parquet"
VCF_1KG = REPO / "data/region_500kb.vcf.gz"
PANEL = Path("/home/yyamada1225/LVexome/demo_delta_mc/thousand_g/"
             "integrated_call_samples_v3.20130502.ALL.panel")

# LIRI paths follow the same convention used in liri_c5_score_test.py and
# friends; override via LIRI_DATA_DIR if your install is elsewhere.
LIRI_BASE = Path(os.environ.get("LIRI_DATA_DIR",
                                 "/mnt/e/LICA-CN-analysis/LIRI-JP_analysis"))
LIRI_VCF = LIRI_BASE / "data/pcawg_germline/LIRI-JP_MHC_germline.vcf.gz"
LIRI_VCF_TO_RK = (LIRI_BASE / "data/pcawg_germline/"
                  "LIRI-JP_rs11509487_genotypes.tsv")
LIRI_CLIN = LIRI_BASE / "data/clinical_data/EGA_clinical_matched.csv"

OUT = REPO / "results"

ANCHOR_POS = 31_366_595          # rs2596542 GRCh37
EAS_POPS = ["CHB", "JPT", "CHS", "CDX", "KHV"]

AXIS_I_COMPS = [4, 6]            # MICA-stable / MICA↓ axis
AXIS_II_COMPS = [5]              # HLA-B / HLA-C-variable, signed-LD reversal


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
def load_1kg_genotypes(vcf: Path, want_pos: set[int]):
    """Phased 1000G VCF -> (samples, positions, haps[len(pos), 2*n_samples])."""
    samples = None; positions = []; haps = []
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            pos = int(f[1])
            ref, alt = f[3], f[4]
            if (len(ref) != 1 or len(alt) != 1
                    or alt.startswith("<") or "," in alt):
                continue
            if pos not in want_pos and pos != ANCHOR_POS:
                continue
            gts = f[9:]
            n = len(gts); hap = np.empty(2 * n, dtype=np.int8)
            ok = True
            for i, g in enumerate(gts):
                if "|" not in g:
                    ok = False; break
                a, b = g.split("|", 1)
                if a == "." or b == ".":
                    hap[2*i] = -1; hap[2*i+1] = -1
                else:
                    hap[2*i] = int(a); hap[2*i+1] = int(b)
            if not ok: continue
            positions.append(pos); haps.append(hap)
    return samples, np.array(positions), np.array(haps, dtype=np.int8)


def load_liri_genotypes(vcf: Path, want_pos: list[int]):
    """Unphased LIRI VCF -> (samples, positions, dosage[len(pos), n_samples])."""
    # bcftools query is fast and produces a regions-restricted dump.
    bed_path = OUT / "_axes_positions.bed"
    bed_path.parent.mkdir(parents=True, exist_ok=True)
    with open(bed_path, "w") as fh:
        for p in want_pos:
            fh.write(f"6\t{p - 1}\t{p}\n")

    samples = subprocess.check_output(
        ["bcftools", "query", "-l", str(vcf)], text=True
    ).strip().split("\n")

    proc = subprocess.run(
        ["bcftools", "query", "-R", str(bed_path),
         "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(vcf)],
        capture_output=True, text=True, timeout=900,
    )
    if proc.returncode != 0:
        print("bcftools STDERR:", proc.stderr[:500])
        raise SystemExit(proc.returncode)

    positions = []
    rows = []
    n_samples = len(samples)
    for line in proc.stdout.strip().split("\n"):
        if not line: continue
        f = line.split("\t")
        pos = int(f[1]); ref, alt = f[2], f[3]
        if "," in alt or alt.startswith("<"): continue
        if len(ref) != 1 or len(alt) != 1: continue
        gts = f[4:]
        if len(gts) != n_samples: continue
        dose = np.full(n_samples, np.nan, dtype=np.float32)
        for i, g in enumerate(gts):
            if "|" in g:    a, b = g.split("|", 1)
            elif "/" in g:  a, b = g.split("/", 1)
            else: continue
            if a == "." or b == ".": continue
            try: dose[i] = int(a) + int(b)
            except ValueError: continue
        positions.append(pos); rows.append(dose)
    return samples, np.array(positions, dtype=int), np.array(rows, dtype=np.float32)


# ---------------------------------------------------------------------------
# Projection
# ---------------------------------------------------------------------------
def project_to_axes(X: np.ndarray, H: np.ndarray) -> np.ndarray:
    """NNLS projection. X is (n_obs, n_snv), H is (k, n_snv).

    Returns W with shape (n_obs, k) such that X[i] ≈ W[i] @ H.
    """
    HT = H.T
    out = np.empty((X.shape[0], H.shape[0]), dtype=np.float32)
    for i in range(X.shape[0]):
        x = X[i]
        valid = ~np.isnan(x)
        if not valid.any():
            out[i] = np.nan; continue
        if valid.all():
            w, _ = nnls(HT, x)
        else:
            w, _ = nnls(HT[valid], x[valid])
        out[i] = w
    return out


def axis_proportions(W: np.ndarray) -> pd.DataFrame:
    """Per-row normalized loadings + axis I / II proportions."""
    K = W.shape[1]
    tot = W.sum(axis=1)
    safe = np.where(tot > 0, tot, np.nan)
    prop = W / safe[:, None]
    df = pd.DataFrame({f"p_c{c}": prop[:, c] for c in range(K)})
    df["axis_I_c4c6"] = prop[:, AXIS_I_COMPS].sum(axis=1)
    df["axis_II_c5"] = prop[:, AXIS_II_COMPS].sum(axis=1)
    df["axis_other"] = 1.0 - df["axis_I_c4c6"] - df["axis_II_c5"]
    df["W_total"] = tot
    return df


# ---------------------------------------------------------------------------
# Stats
# ---------------------------------------------------------------------------
def mwu(a: np.ndarray, b: np.ndarray):
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 3 or len(b) < 3:
        return float("nan"), float("nan")
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    return float(u), float(p)


def desc(series: np.ndarray) -> dict:
    s = series[~np.isnan(series)]
    if len(s) == 0:
        return {"n": 0, "mean": float("nan"), "median": float("nan"),
                "q25": float("nan"), "q75": float("nan")}
    return {
        "n": int(len(s)),
        "mean": float(np.mean(s)),
        "median": float(np.median(s)),
        "q25": float(np.quantile(s, 0.25)),
        "q75": float(np.quantile(s, 0.75)),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    sys.stdout.reconfigure(line_buffering=True)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # ---- NMF H ----
    H_df = pd.read_parquet(H_PARQUET)
    nmf_pos = np.array([int(c) for c in H_df.columns], dtype=int)
    H = H_df.values.astype(np.float32)                       # (k, n_snv)
    K = H.shape[0]
    nmf_pos_set = set(int(p) for p in nmf_pos)
    print(f"NMF H: k={K}, n_SNV={H.shape[1]}, "
          f"window {nmf_pos.min():,}-{nmf_pos.max():,}")

    # ---- 1000G phased haps ----
    print("\n[1000G] loading VCF (restricted to NMF positions + anchor)...")
    samples_1k, pos_1k, haps_1k = load_1kg_genotypes(VCF_1KG, nmf_pos_set)
    print(f"  samples={len(samples_1k)}, positions in VCF={len(pos_1k)}, "
          f"haps={haps_1k.shape[1]}  ({time.time()-t0:.0f}s)")

    panel = pd.read_csv(PANEL, sep="\t")
    pop_of = dict(zip(panel["sample"], panel["pop"]))

    # Anchor row + position alignment
    if ANCHOR_POS not in pos_1k:
        raise SystemExit("anchor not found in 1000G VCF")
    anchor_idx_1k = int(np.where(pos_1k == ANCHOR_POS)[0][0])
    anchor_alleles_1k = haps_1k[anchor_idx_1k]               # (2*n_samples,)

    # Per-NMF-position 1000G dosage / hap matrix (drop anchor row, align)
    pos_to_row_1k = {int(p): i for i, p in enumerate(pos_1k)}
    keep_idx = [pos_to_row_1k[p] for p in nmf_pos if p in pos_to_row_1k]
    used_nmf_mask = np.array([p in pos_to_row_1k for p in nmf_pos], dtype=bool)
    H_used = H[:, used_nmf_mask]                              # (k, n_used)
    used_pos = nmf_pos[used_nmf_mask]
    haps_1k_aligned = haps_1k[keep_idx]                       # (n_used, 2*n_samples)
    print(f"  NMF SNVs aligned to 1000G: "
          f"{used_nmf_mask.sum()}/{len(nmf_pos)} ({time.time()-t0:.0f}s)")

    # Per-hap 1000G EAS  --------------------------------------------------
    pop_arr = np.array([pop_of.get(s, "UNK") for s in samples_1k])
    is_eas_sample = np.isin(pop_arr, EAS_POPS)
    eas_sample_idx = np.where(is_eas_sample)[0]
    # Per-haplotype: 2 haps per sample
    hap_sample_idx = np.repeat(np.arange(len(samples_1k)), 2)
    hap_pop = np.repeat(pop_arr, 2)
    hap_eas_mask = np.isin(hap_pop, EAS_POPS)
    eas_hap_idx = np.where(hap_eas_mask)[0]

    # Build per-hap dosage matrix (haps × SNV)
    X_hap_eas = haps_1k_aligned[:, eas_hap_idx].T.astype(np.float32)
    # Missing genotypes (-1) -> NaN so projection skips them.
    X_hap_eas = np.where(X_hap_eas < 0, np.nan, X_hap_eas)
    is_T_carrier_hap_eas = (anchor_alleles_1k[eas_hap_idx] == 1)
    print(f"  1000G EAS haps: total={len(eas_hap_idx)}, "
          f"T-carrier={is_T_carrier_hap_eas.sum()}, "
          f"non-carrier={(~is_T_carrier_hap_eas).sum()}")

    # Per-sample diploid dosage (EAS)
    n_eas_samples = len(eas_sample_idx)
    X_dip_eas = np.empty((n_eas_samples, H_used.shape[1]), dtype=np.float32)
    anchor_dose_eas = np.zeros(n_eas_samples, dtype=int)
    for k, si in enumerate(eas_sample_idx):
        h0 = haps_1k_aligned[:, 2*si].astype(np.float32)
        h1 = haps_1k_aligned[:, 2*si+1].astype(np.float32)
        h0 = np.where(h0 < 0, np.nan, h0)
        h1 = np.where(h1 < 0, np.nan, h1)
        X_dip_eas[k] = np.nansum(np.stack([h0, h1]), axis=0)
        # Mark fully-missing as NaN
        both_missing = np.isnan(h0) & np.isnan(h1)
        X_dip_eas[k, both_missing] = np.nan
        a0 = anchor_alleles_1k[2*si]; a1 = anchor_alleles_1k[2*si+1]
        if a0 >= 0 and a1 >= 0:
            anchor_dose_eas[k] = int(a0) + int(a1)
        else:
            anchor_dose_eas[k] = -1

    # Project per-hap and per-sample
    print(f"\n[1000G] projecting per-hap onto H ({time.time()-t0:.0f}s)...")
    W_hap_eas = project_to_axes(X_hap_eas, H_used)
    print(f"[1000G] projecting per-sample onto 2H ({time.time()-t0:.0f}s)...")
    W_dip_eas = project_to_axes(X_dip_eas, 2 * H_used)

    df_hap_eas = axis_proportions(W_hap_eas)
    df_hap_eas["cohort"] = "1000G_EAS"
    df_hap_eas["pop"] = hap_pop[eas_hap_idx]
    df_hap_eas["anchor_carrier"] = np.where(is_T_carrier_hap_eas,
                                              "T_carrier", "non_carrier")

    df_dip_eas = axis_proportions(W_dip_eas)
    df_dip_eas["cohort"] = "1000G_EAS"
    df_dip_eas["pop"] = pop_arr[eas_sample_idx]
    df_dip_eas["anchor_dose"] = anchor_dose_eas
    df_dip_eas["anchor_carrier"] = np.where(
        anchor_dose_eas <= 0, "non_carrier",
        np.where(anchor_dose_eas == 1, "het", "hom_T"))

    # ---- LIRI-JP ----
    print(f"\n[LIRI] loading VCF ({time.time()-t0:.0f}s)...")
    if not LIRI_VCF.exists():
        print(f"  WARN: LIRI VCF not found at {LIRI_VCF} -- skipping LIRI side.")
        df_hap_liri = pd.DataFrame()
        df_dip_liri = pd.DataFrame()
    else:
        want_positions = sorted(set(int(p) for p in used_pos) | {ANCHOR_POS})
        samp_liri, pos_liri, dose_liri = load_liri_genotypes(
            LIRI_VCF, want_positions)
        print(f"  LIRI: samples={len(samp_liri)}, positions returned={len(pos_liri)}")
        # Align LIRI positions to used_pos
        liri_pos_to_row = {int(p): i for i, p in enumerate(pos_liri)}
        aligned_rows = [liri_pos_to_row.get(int(p), -1) for p in used_pos]
        n_used = H_used.shape[1]
        X_dip_liri = np.full((len(samp_liri), n_used), np.nan, dtype=np.float32)
        for j, ri in enumerate(aligned_rows):
            if ri >= 0:
                X_dip_liri[:, j] = dose_liri[ri]
        # Anchor dose per LIRI sample
        if ANCHOR_POS in liri_pos_to_row:
            anchor_dose_liri = dose_liri[liri_pos_to_row[ANCHOR_POS]].copy()
        else:
            print("  WARN: rs2596542 not in LIRI VCF; treating all as carriers.")
            anchor_dose_liri = np.full(len(samp_liri), np.nan)

        # Drop samples with no usable data
        usable = (~np.isnan(X_dip_liri)).sum(axis=1) >= max(1, int(0.5 * n_used))
        print(f"  LIRI samples with >=50% usable SNVs: {usable.sum()}/{usable.size}")

        print(f"[LIRI] projecting per-sample onto 2H ({time.time()-t0:.0f}s)...")
        W_dip_liri_full = project_to_axes(X_dip_liri, 2 * H_used)
        df_dip_liri = axis_proportions(W_dip_liri_full)
        df_dip_liri["sample"] = samp_liri
        df_dip_liri["cohort"] = "LIRI"
        df_dip_liri["anchor_dose"] = [
            int(d) if not np.isnan(d) else -1 for d in anchor_dose_liri
        ]
        df_dip_liri["anchor_carrier"] = np.where(
            df_dip_liri["anchor_dose"].values <= 0, "non_carrier",
            np.where(df_dip_liri["anchor_dose"].values == 1, "het", "hom_T"))
        df_dip_liri["usable_pct"] = (~np.isnan(X_dip_liri)).mean(axis=1)
        df_dip_liri = df_dip_liri[usable].reset_index(drop=True)

        # ---- Per-hap proxy for LIRI: split dosage in half across two pseudo-haps.
        # This is a heuristic for unphased data; mean of two haps == diploid/2.
        # Both pseudo-haps inherit the projection of dosage/2, so per-hap
        # proportions equal per-sample proportions but doubled in count.
        # We instead emit a per-hap table that EXPANDS each diploid into two
        # rows with anchor_carrier set by the diploid genotype (het -> 1 T-hap,
        # 1 non-T-hap; hom -> 2 T-haps).
        rows = []
        for idx, row in df_dip_liri.iterrows():
            dose = row["anchor_dose"]
            for k_hap in range(2):
                if dose == 2:    carrier = "T_carrier"
                elif dose == 0:  carrier = "non_carrier"
                elif dose == 1:  carrier = "T_carrier" if k_hap == 0 else "non_carrier"
                else:            carrier = "unknown"
                rows.append({
                    **{c: row[c] for c in df_dip_liri.columns
                       if c not in ("anchor_dose", "anchor_carrier",
                                    "usable_pct", "sample", "cohort")},
                    "cohort": "LIRI", "pop": "JPT_HCC",
                    "anchor_carrier": carrier,
                })
        df_hap_liri = pd.DataFrame(rows)

    # ---- Save per-hap and per-sample tables ----
    if not df_hap_liri.empty:
        df_hap = pd.concat([df_hap_eas, df_hap_liri], ignore_index=True)
    else:
        df_hap = df_hap_eas
    if not df_dip_liri.empty:
        df_dip = pd.concat([df_dip_eas, df_dip_liri], ignore_index=True)
    else:
        df_dip = df_dip_eas
    df_hap.to_csv(OUT / "cohort_axis_proportions_perhap.csv", index=False)
    df_dip.to_csv(OUT / "cohort_axis_proportions_persample.csv", index=False)
    print(f"\nwrote {OUT / 'cohort_axis_proportions_perhap.csv'}")
    print(f"wrote {OUT / 'cohort_axis_proportions_persample.csv'}")

    # ---- Optional: LIRI viral_status merge ----
    df_dip_viral = df_dip_liri.copy() if not df_dip_liri.empty else pd.DataFrame()
    if LIRI_VCF_TO_RK.exists() and LIRI_CLIN.exists() and not df_dip_viral.empty:
        umap = pd.read_csv(LIRI_VCF_TO_RK, sep="\t",
                            usecols=["liriju_sample_id", "vcf_sample_id"])
        umap = umap.rename(columns={"vcf_sample_id": "sample",
                                     "liriju_sample_id": "rk_id"})
        clin = pd.read_csv(LIRI_CLIN).rename(
            columns={"clinical_sample_id": "rk_id"})
        df_dip_viral = df_dip_viral.merge(umap, on="sample", how="left") \
                                    .merge(clin[["rk_id", "viral_status"]],
                                           on="rk_id", how="left")
        df_dip_viral.to_csv(
            OUT / "cohort_axis_proportions_persample_with_viral.csv",
            index=False)
        print(f"wrote {OUT / 'cohort_axis_proportions_persample_with_viral.csv'}")

    # ---- Statistics ----
    stats = []

    def add_stat(scope, axis, group_a, a_vals, group_b, b_vals):
        u, p = mwu(np.asarray(a_vals, dtype=float),
                   np.asarray(b_vals, dtype=float))
        da = desc(np.asarray(a_vals, dtype=float))
        db = desc(np.asarray(b_vals, dtype=float))
        stats.append({
            "scope": scope, "axis": axis,
            "group_A": group_a, "n_A": da["n"],
            "mean_A": da["mean"], "median_A": da["median"],
            "group_B": group_b, "n_B": db["n"],
            "mean_B": db["mean"], "median_B": db["median"],
            "diff_mean_A_minus_B": da["mean"] - db["mean"],
            "MWU_U": u, "MWU_p": p,
        })

    for axis_name, axis_col in [("Axis_I_c4c6", "axis_I_c4c6"),
                                  ("Axis_II_c5", "axis_II_c5")]:
        # per-hap T-carrier comparison: 1000G EAS T vs LIRI T
        a = df_hap_eas.loc[df_hap_eas["anchor_carrier"] == "T_carrier",
                            axis_col].values
        if not df_hap_liri.empty:
            b = df_hap_liri.loc[df_hap_liri["anchor_carrier"] == "T_carrier",
                                 axis_col].values
            add_stat("perhap_Tcarrier", axis_name,
                     "1000G_EAS_T", a, "LIRI_T", b)
            # per-hap non-carrier
            a2 = df_hap_eas.loc[df_hap_eas["anchor_carrier"] == "non_carrier",
                                 axis_col].values
            b2 = df_hap_liri.loc[df_hap_liri["anchor_carrier"] == "non_carrier",
                                  axis_col].values
            add_stat("perhap_nonCarrier", axis_name,
                     "1000G_EAS_nonT", a2, "LIRI_nonT", b2)
        # per-hap, within each cohort: carrier vs non-carrier
        a3 = df_hap_eas.loc[df_hap_eas["anchor_carrier"] == "T_carrier",
                             axis_col].values
        b3 = df_hap_eas.loc[df_hap_eas["anchor_carrier"] == "non_carrier",
                             axis_col].values
        add_stat("perhap_carrierVsNon_1000G_EAS", axis_name,
                 "1000G_EAS_T", a3, "1000G_EAS_nonT", b3)
        if not df_hap_liri.empty:
            a4 = df_hap_liri.loc[df_hap_liri["anchor_carrier"] == "T_carrier",
                                  axis_col].values
            b4 = df_hap_liri.loc[df_hap_liri["anchor_carrier"] == "non_carrier",
                                  axis_col].values
            add_stat("perhap_carrierVsNon_LIRI", axis_name,
                     "LIRI_T", a4, "LIRI_nonT", b4)
        # per-sample: 1000G EAS carriers (dose>=1) vs LIRI carriers
        a5 = df_dip_eas.loc[df_dip_eas["anchor_dose"] >= 1, axis_col].values
        b5 = (df_dip_liri.loc[df_dip_liri["anchor_dose"] >= 1, axis_col].values
               if not df_dip_liri.empty else np.array([]))
        if not df_dip_liri.empty:
            add_stat("persample_carrier", axis_name,
                     "1000G_EAS_carrier", a5, "LIRI_carrier", b5)
            a6 = df_dip_eas.loc[df_dip_eas["anchor_dose"] == 0, axis_col].values
            b6 = df_dip_liri.loc[df_dip_liri["anchor_dose"] == 0,
                                  axis_col].values
            add_stat("persample_nonCarrier", axis_name,
                     "1000G_EAS_nonCarrier", a6, "LIRI_nonCarrier", b6)
        # JPT-specific
        a7 = df_dip_eas.loc[(df_dip_eas["pop"] == "JPT")
                              & (df_dip_eas["anchor_dose"] >= 1),
                              axis_col].values
        if not df_dip_liri.empty:
            add_stat("persample_carrier_JPTvsLIRI", axis_name,
                     "1000G_JPT_carrier", a7,
                     "LIRI_carrier",
                     df_dip_liri.loc[df_dip_liri["anchor_dose"] >= 1,
                                       axis_col].values)
        # LIRI viral_status stratification (carriers only)
        if not df_dip_viral.empty and "viral_status" in df_dip_viral.columns:
            car = df_dip_viral[df_dip_viral["anchor_dose"] >= 1]
            for vs_a, vs_b in [("HBV", "HCV"), ("HBV", "NBNC"), ("HCV", "NBNC")]:
                aa = car.loc[car["viral_status"] == vs_a, axis_col].values
                bb = car.loc[car["viral_status"] == vs_b, axis_col].values
                add_stat(f"LIRI_carrier_{vs_a}_vs_{vs_b}",
                         axis_name, vs_a, aa, vs_b, bb)

    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(OUT / "cohort_axis_proportions_stats.csv", index=False)
    print(f"wrote {OUT / 'cohort_axis_proportions_stats.csv'}")
    print("\n=== Mann-Whitney U comparisons ===")
    print(stats_df.to_string(index=False, formatters={
        "mean_A": "{:.3f}".format, "median_A": "{:.3f}".format,
        "mean_B": "{:.3f}".format, "median_B": "{:.3f}".format,
        "diff_mean_A_minus_B": "{:+.3f}".format,
        "MWU_U": "{:.0f}".format, "MWU_p": "{:.3g}".format,
    }))

    # ---- Figure ----
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
    })

    def violin(ax, groups: list[tuple[str, np.ndarray, str]], ylabel: str,
               title: str):
        positions = list(range(len(groups)))
        data = [g[1][~np.isnan(g[1])] for g in groups]
        # Guard against empty groups
        plot_data = [d if len(d) else np.array([np.nan]) for d in data]
        parts = ax.violinplot(plot_data, positions=positions,
                                showmeans=True, showmedians=False, widths=0.78)
        for body, (_, _, col) in zip(parts["bodies"], groups):
            body.set_facecolor(col); body.set_alpha(0.45); body.set_edgecolor("#222")
        if "cmeans" in parts:
            parts["cmeans"].set_color("#000")
        rng = np.random.default_rng(0)
        for i, (lab, vals, col) in enumerate(groups):
            v = vals[~np.isnan(vals)]
            if len(v) == 0: continue
            xj = i + rng.uniform(-0.12, 0.12, len(v))
            ax.scatter(xj, v, s=6, color=col, alpha=0.3, edgecolor="none")
        labels = [f"{g[0]}\nn={int(np.sum(~np.isnan(g[1])))}" for g in groups]
        ax.set_xticks(positions); ax.set_xticklabels(labels, fontsize=7.5)
        ax.set_ylabel(ylabel); ax.set_title(title, loc="left", fontsize=10.5)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(axis="y", linestyle=":", alpha=0.3)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # Per-haplotype, Axis I
    grp_hap_I = [
        ("1000G EAS\nT-carrier",
         df_hap_eas.loc[df_hap_eas["anchor_carrier"] == "T_carrier",
                         "axis_I_c4c6"].values, "#1f4e8c"),
        ("1000G EAS\nnon-carrier",
         df_hap_eas.loc[df_hap_eas["anchor_carrier"] == "non_carrier",
                         "axis_I_c4c6"].values, "#9ab8d8"),
    ]
    if not df_hap_liri.empty:
        grp_hap_I += [
            ("LIRI HCC\nT-carrier (proxy)",
             df_hap_liri.loc[df_hap_liri["anchor_carrier"] == "T_carrier",
                              "axis_I_c4c6"].values, "#d65454"),
            ("LIRI HCC\nnon-carrier (proxy)",
             df_hap_liri.loc[df_hap_liri["anchor_carrier"] == "non_carrier",
                              "axis_I_c4c6"].values, "#e9a8a8"),
        ]
    violin(axes[0, 0], grp_hap_I, "Axis I proportion (c4 + c6)",
           "(a)  per-haplotype  —  Axis I  (MICA-stable / MICA↓)")

    grp_hap_II = [
        ("1000G EAS\nT-carrier",
         df_hap_eas.loc[df_hap_eas["anchor_carrier"] == "T_carrier",
                         "axis_II_c5"].values, "#1f4e8c"),
        ("1000G EAS\nnon-carrier",
         df_hap_eas.loc[df_hap_eas["anchor_carrier"] == "non_carrier",
                         "axis_II_c5"].values, "#9ab8d8"),
    ]
    if not df_hap_liri.empty:
        grp_hap_II += [
            ("LIRI HCC\nT-carrier (proxy)",
             df_hap_liri.loc[df_hap_liri["anchor_carrier"] == "T_carrier",
                              "axis_II_c5"].values, "#d65454"),
            ("LIRI HCC\nnon-carrier (proxy)",
             df_hap_liri.loc[df_hap_liri["anchor_carrier"] == "non_carrier",
                              "axis_II_c5"].values, "#e9a8a8"),
        ]
    violin(axes[0, 1], grp_hap_II, "Axis II proportion (c5)",
           "(b)  per-haplotype  —  Axis II  (HLA-B/HLA-C-variable; signed-LD reversal)")

    # Per-sample, Axis I (anchor dose 0 vs 1 vs 2)
    def by_dose(df, axis_col, doses, color_map):
        groups = []
        for d in doses:
            sel = (df["anchor_dose"] == d).values
            groups.append((f"dose {d}", df.loc[sel, axis_col].values,
                            color_map[d]))
        return groups

    color_eas = {0: "#9ab8d8", 1: "#3070b8", 2: "#1f4e8c"}
    color_liri = {0: "#e9a8a8", 1: "#d65454", 2: "#a23030"}

    grp_dip_I = []
    grp_dip_II = []
    for d, lab, col_eas in [(0, "EAS dose 0", color_eas[0]),
                              (1, "EAS dose 1", color_eas[1]),
                              (2, "EAS dose 2", color_eas[2])]:
        sel = (df_dip_eas["anchor_dose"] == d).values
        grp_dip_I.append((lab,
                          df_dip_eas.loc[sel, "axis_I_c4c6"].values, col_eas))
        grp_dip_II.append((lab,
                           df_dip_eas.loc[sel, "axis_II_c5"].values, col_eas))
    if not df_dip_liri.empty:
        for d, lab, col_l in [(0, "LIRI dose 0", color_liri[0]),
                                (1, "LIRI dose 1", color_liri[1]),
                                (2, "LIRI dose 2", color_liri[2])]:
            sel = (df_dip_liri["anchor_dose"] == d).values
            grp_dip_I.append((lab,
                              df_dip_liri.loc[sel, "axis_I_c4c6"].values,
                              col_l))
            grp_dip_II.append((lab,
                               df_dip_liri.loc[sel, "axis_II_c5"].values,
                               col_l))

    violin(axes[1, 0], grp_dip_I, "Axis I proportion (c4 + c6)",
           "(c)  per-sample  —  Axis I by rs2596542-T dose")
    violin(axes[1, 1], grp_dip_II, "Axis II proportion (c5)",
           "(d)  per-sample  —  Axis II by rs2596542-T dose")

    fig.suptitle(
        "NMF axis proportions: 1000G EAS vs LIRI-JP HCC  —  "
        "carrier / non-carrier breakdown\n"
        "Axis I = (c4 + c6) / Σ_k Wₖ  (MICA-stable),  "
        "Axis II = c5 / Σ_k Wₖ  (HLA-B/HLA-C-variable)",
        fontsize=11.5, y=1.01,
    )
    fig.tight_layout()
    out_fig = OUT / "figure_cohort_axis_proportions"
    fig.savefig(out_fig.with_suffix(".pdf"))
    fig.savefig(out_fig.with_suffix(".png"), dpi=300)
    print(f"wrote {out_fig}.pdf / .png")
    plt.close(fig)

    print(f"\nTotal: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
