"""MAF/distance-matched permutation test for c5 enrichment.

For each c5 top-5% SNV, draw a "matched control" from the universe with:
  - similar MAF (within carriers, +/- 1 decile)
  - similar distance to anchor (+/- 1 decile)

Repeat 1000 permutations. Compare observed fold enrichment to null
distribution for:
  - HLA-B Whole_Blood eQTL set
  - HLA-C Whole_Blood eQTL set
  - MICB Whole_Blood eQTL set
  - any-pair C-flip set
  - any-pair r-flip set

Outputs:
  results/c5_permutation_null.csv
  results/figure_c5_permutation.{pdf,png}
"""

from __future__ import annotations

import json
import sys
import time
import urllib.request
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
H_PARQUET = REPO / "results/nmf_H_26_k8.parquet"
EQTL_CSV = REPO / "results/eqtl_targets_per_gene.csv"
REGION_CSV = REPO / "results/region_per_partner_26.csv"
OUT = REPO / "results"

ANCHOR_POS = 31_366_595
TARGET_COMPONENT = 5
N_PERMUTE = 1000
NMF_TOP_FRAC = 0.05
RNG = np.random.default_rng(0)
SUPER_POPS = ["AFR", "EUR", "EAS", "SAS", "AMR"]


def main():
    sys.stdout.reconfigure(line_buffering=True)
    H = pd.read_parquet(H_PARQUET)
    snv_pos = np.array([int(c) for c in H.columns])
    h_mat = H.values
    n = len(snv_pos)
    top_n = max(1, int(NMF_TOP_FRAC * n))
    print(f"NMF: K={h_mat.shape[0]}, n SNVs={n}, top-5% size={top_n}")

    # Need MAF per SNV in carrier set. Re-derive (or use H column heuristic).
    # We didn't save MAF earlier — recompute by reading region CSV (median across pops)
    # or just use H column means. Actually the cleanest is to reload region_per_pop
    # and average p_b across pops.
    print("loading region_per_pop_26.parquet for MAF...")
    rpp = pd.read_parquet(REPO / "results/region_per_pop_26.parquet")
    # Use EAS_pooled p_b as a proxy MAF (or all pops mean)
    # Actually, p_b is per-pop, varies. Use mean across all 5 super_pop pooled.
    pooled = rpp[rpp["pop"].isin(["AFR_pooled","EUR_pooled","EAS_pooled","SAS_pooled","AMR_pooled"])]
    maf_per_pos = pooled.groupby("partner_pos")["p_b"].mean()
    # Map back to snv_pos
    maf_arr = np.array([maf_per_pos.get(int(p), float("nan")) for p in snv_pos])
    valid = ~np.isnan(maf_arr)
    print(f"  MAF assigned for {valid.sum()}/{n}")

    # Distance from anchor
    dist_arr = np.abs(snv_pos - ANCHOR_POS).astype(float)

    # Stratify by MAF decile × distance decile (10×10 = 100 strata)
    maf_dec = pd.qcut(maf_arr[valid], 10, labels=False, duplicates="drop")
    dist_dec = pd.qcut(dist_arr[valid], 10, labels=False, duplicates="drop")
    stratum_full = np.full(n, -1, dtype=int)
    stratum_full[valid] = maf_dec * 10 + dist_dec
    print(f"  unique strata: {len(set(stratum_full[stratum_full>=0]))}")

    # c5 top-5%
    c5_top_idx = np.argsort(h_mat[TARGET_COMPONENT])[-top_n:]
    c5_top_pos = snv_pos[c5_top_idx]
    c5_top_set = set(int(p) for p in c5_top_pos)
    c5_strata = stratum_full[c5_top_idx]

    # Build target-gene eQTL sets (in window, b37) from saved CSV
    eqtl_df = pd.read_csv(EQTL_CSV)
    eqtl_sets = {
        gene: set(eqtl_df[eqtl_df["gene"] == gene]["b37_pos"].astype(int))
        for gene in eqtl_df["gene"].unique()
    }
    universe = set(int(p) for p in snv_pos)
    for g in eqtl_sets:
        eqtl_sets[g] &= universe

    # Build flip sets from region CSV
    region = pd.read_csv(REGION_CSV)
    cflip_set = set(region.loc[region["any_C_flip"] == True, "partner_pos"].astype(int)) & universe
    rmat = region[[f"median_r_{s}" for s in SUPER_POPS]].astype(float).values
    nrf = np.zeros(len(region), dtype=int)
    np_pairs = np.zeros(len(region), dtype=int)
    for i, j in combinations(range(len(SUPER_POPS)), 2):
        ri = rmat[:, i]; rj = rmat[:, j]
        v = ~(np.isnan(ri) | np.isnan(rj))
        np_pairs += v.astype(int)
        flip = v & ((np.sign(ri) * np.sign(rj)) < 0)
        nrf += flip.astype(int)
    rflip_set = set(region.loc[nrf >= 1, "partner_pos"].astype(int)) & universe

    # Test sets
    test_sets = {
        "HLA-B": eqtl_sets.get("HLA-B", set()),
        "HLA-C": eqtl_sets.get("HLA-C", set()),
        "MICB":  eqtl_sets.get("MICB",  set()),
        "MICA":  eqtl_sets.get("MICA",  set()),
        "C_flip": cflip_set,
        "r_flip": rflip_set,
    }
    print(f"\ntest set sizes:")
    for k, s in test_sets.items():
        print(f"  {k}: {len(s)}")

    # Observed enrichment
    print(f"\nobserved c5 top-{top_n} overlap:")
    obs = {}
    for k, s in test_sets.items():
        x = len(c5_top_set & s)
        exp = top_n * len(s) / len(universe)
        obs[k] = {"x": x, "exp": exp, "fold": x / exp if exp > 0 else float("nan")}
        print(f"  {k}: x={x}, expected={exp:.1f}, fold={obs[k]['fold']:.2f}")

    # Build stratum index for permutation sampling
    stratum_to_indices = {}
    for j in np.where(valid)[0]:
        stratum_to_indices.setdefault(int(stratum_full[j]), []).append(j)
    for s in stratum_to_indices:
        stratum_to_indices[s] = np.array(stratum_to_indices[s])

    # Permutation
    print(f"\nrunning {N_PERMUTE} matched permutations...")
    t0 = time.time()
    null_folds = {k: np.zeros(N_PERMUTE) for k in test_sets}
    for p in range(N_PERMUTE):
        sampled_idx = []
        for s in c5_strata:
            if s < 0:
                # Random fallback
                sampled_idx.append(RNG.integers(0, n))
            else:
                pool = stratum_to_indices.get(int(s), None)
                if pool is None or len(pool) == 0:
                    sampled_idx.append(RNG.integers(0, n))
                else:
                    sampled_idx.append(RNG.choice(pool))
        sampled_pos = set(int(snv_pos[i]) for i in sampled_idx)
        for k, s in test_sets.items():
            x = len(sampled_pos & s)
            exp = len(sampled_pos) * len(s) / len(universe)
            null_folds[k][p] = x / exp if exp > 0 else 0.0
        if (p + 1) % 200 == 0:
            print(f"  ...{p+1}/{N_PERMUTE} ({time.time()-t0:.0f}s)")

    # Permutation p-value (one-sided)
    print(f"\nObserved fold vs MAF×distance-matched null (1000 perms):")
    rows = []
    for k in test_sets:
        obs_fold = obs[k]["fold"]
        null = null_folds[k]
        p_perm = float((null >= obs_fold).sum() / N_PERMUTE)
        z = (obs_fold - null.mean()) / (null.std() if null.std() > 0 else 1.0)
        rows.append({
            "test_set": k,
            "obs_fold": obs_fold,
            "null_mean": float(null.mean()),
            "null_sd": float(null.std()),
            "z_score": z,
            "perm_p_one_sided": p_perm,
            "n_obs_overlap": obs[k]["x"],
        })
        print(f"  {k:8s}: obs={obs_fold:.2f}  null mean={null.mean():.2f}±{null.std():.2f}  "
              f"z={z:+.2f}  perm p={p_perm:.4f}")
    pd.DataFrame(rows).to_csv(OUT / "c5_permutation_null.csv", index=False)

    # ===== Figure =====
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.5))
    test_order = ["HLA-B", "HLA-C", "MICB", "MICA", "C_flip", "r_flip"]
    palette = {"HLA-B":"#5a4fcf","HLA-C":"#3a8c5a","MICB":"#b5722d",
                "MICA":"#d65454","C_flip":"#7e3da7","r_flip":"#1f4e8c"}
    for ax, k in zip(axes.flat, test_order):
        null = null_folds[k]
        obs_fold = obs[k]["fold"]
        ax.hist(null, bins=40, color=palette[k], alpha=0.55, edgecolor="#222",
                 linewidth=0.4)
        ax.axvline(obs_fold, color="#d65454", lw=1.6, ls="-",
                    label=f"obs={obs_fold:.2f}")
        ax.axvline(null.mean(), color="#222", lw=0.7, ls="--",
                    label=f"null mean={null.mean():.2f}")
        p_p = (null >= obs_fold).sum() / N_PERMUTE
        ax.set_title(f"{k}  (perm p = {p_p:.4f})", loc="left", fontsize=10)
        ax.set_xlabel("fold enrichment")
        ax.set_ylabel("permutations")
        ax.legend(fontsize=7.5, loc="upper right")
    fig.suptitle(
        "MAF × distance-matched permutation null for c5 top-5% SNV enrichment\n"
        "(HLA-B-variable axis, 1000 perms, anchor rs2596542-T)",
        fontsize=11, y=1.00,
    )
    fig.tight_layout()
    out = OUT / "figure_c5_permutation"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
