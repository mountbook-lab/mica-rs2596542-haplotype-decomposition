"""Per-component MAF×distance matched permutation (k=8, all 8 components).

For each component:
  observed = overlap of top-5% SNVs with each test set
  permutation: 10,000 matched draws (same MAF bin, same distance bin)
  empirical p (one-sided), z-score, fold vs matched

Bins:
  MAF (overall p_b mean across pooled super-pops):
    [0, 0.01), [0.01, 0.05), [0.05, 0.10), [0.10, 0.20),
    [0.20, 0.30), [0.30, 1.0]
  distance from rs2596542 (bp):
    [0, 25k), [25k, 50k), [50k, 100k), [100k, 150k),
    [150k, 200k), [200k, 250k]

Test sets:
  r-flip, C-flip, class-flip, MICA, HLA-B, HLA-C, MICB

Output:
  results/permutation_all_components.csv
  results/figure_permutation_all_components.{pdf,png}
"""

from __future__ import annotations

import sys
import time
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
H_PARQUET = REPO / "results/nmf_H_26_k8.parquet"
EQTL_CSV = REPO / "results/eqtl_targets_per_gene.csv"
REGION_CSV = REPO / "results/region_per_partner_26.csv"
RPP_PARQUET = REPO / "results/region_per_pop_26.parquet"
OUT = REPO / "results"

ANCHOR_POS = 31_366_595
SUPER_POPS = ["AFR", "EUR", "EAS", "SAS", "AMR"]
NMF_TOP_FRAC = 0.05
N_PERMUTE = 10_000
RNG = np.random.default_rng(0)

MAF_BINS = [0.0, 0.01, 0.05, 0.10, 0.20, 0.30, 1.01]
DIST_BINS = [0, 25_000, 50_000, 100_000, 150_000, 200_000, 260_000]


def main():
    sys.stdout.reconfigure(line_buffering=True)
    H = pd.read_parquet(H_PARQUET)
    snv_pos = np.array([int(c) for c in H.columns])
    h_mat = H.values
    K = h_mat.shape[0]
    n = len(snv_pos)
    top_n = max(1, int(NMF_TOP_FRAC * n))

    # Build MAF for each SNV using pooled super-pop p_b mean
    print("loading MAF + distance ...")
    rpp = pd.read_parquet(RPP_PARQUET)
    pooled = rpp[rpp["pop"].isin([f"{s}_pooled" for s in SUPER_POPS])]
    maf_per_pos = pooled.groupby("partner_pos")["p_b"].mean()
    maf_arr = np.array([maf_per_pos.get(int(p), float("nan")) for p in snv_pos])
    dist_arr = np.abs(snv_pos - ANCHOR_POS)

    # Bin assignment
    maf_bin = np.digitize(maf_arr, MAF_BINS) - 1
    dist_bin = np.digitize(dist_arr, DIST_BINS) - 1
    stratum = maf_bin * 100 + dist_bin
    valid = (~np.isnan(maf_arr)) & (maf_bin >= 0) & (dist_bin >= 0)
    print(f"  valid SNVs: {valid.sum()}/{n}")

    stratum_to_idx = {}
    for j in np.where(valid)[0]:
        stratum_to_idx.setdefault(int(stratum[j]), []).append(j)
    for s in stratum_to_idx:
        stratum_to_idx[s] = np.array(stratum_to_idx[s])

    # Test sets
    universe = set(int(p) for p in snv_pos)
    eqtl = pd.read_csv(EQTL_CSV)
    eqtl_sets = {g: set(eqtl[eqtl["gene"] == g]["b37_pos"].astype(int)) & universe
                 for g in eqtl["gene"].unique()}
    region = pd.read_csv(REGION_CSV)
    cflip_set = set(region.loc[region["any_C_flip"] == True,
                                  "partner_pos"].astype(int)) & universe
    rmat = region[[f"median_r_{s}" for s in SUPER_POPS]].astype(float).values
    nrf = np.zeros(len(region), dtype=int)
    for i, j in combinations(range(len(SUPER_POPS)), 2):
        ri = rmat[:, i]; rj = rmat[:, j]
        v = ~(np.isnan(ri) | np.isnan(rj))
        nrf += (v & ((np.sign(ri) * np.sign(rj)) < 0)).astype(int)
    rflip_set = set(region.loc[nrf >= 1, "partner_pos"].astype(int)) & universe
    classflip_set = set(region.loc[region["any_class_flip"] == True,
                                      "partner_pos"].astype(int)) & universe

    test_sets = {
        "r_flip": rflip_set, "C_flip": cflip_set, "class_flip": classflip_set,
        "MICA": eqtl_sets.get("MICA", set()),
        "MICB": eqtl_sets.get("MICB", set()),
        "HLA-B": eqtl_sets.get("HLA-B", set()),
        "HLA-C": eqtl_sets.get("HLA-C", set()),
    }

    # Per-component permutation
    rows = []
    for c in range(K):
        h = h_mat[c]
        top_idx = np.argsort(h)[-top_n:]
        c_strata = stratum[top_idx]
        c_top_set = set(int(p) for p in snv_pos[top_idx])
        # Observed
        obs_x = {k: len(c_top_set & ts) for k, ts in test_sets.items()}

        # Build per-stratum pool
        n_perm_eff = N_PERMUTE
        null_x = {k: np.zeros(n_perm_eff, dtype=int) for k in test_sets}

        # For each permutation, sample one matched SNV per top SNV from same stratum
        t0 = time.time()
        for p in range(n_perm_eff):
            sampled_idx = []
            for s in c_strata:
                pool = stratum_to_idx.get(int(s), None)
                if pool is None or len(pool) == 0:
                    sampled_idx.append(RNG.integers(0, n))
                else:
                    sampled_idx.append(RNG.choice(pool))
            sampled_pos = set(int(snv_pos[i]) for i in sampled_idx)
            for k, ts in test_sets.items():
                null_x[k][p] = len(sampled_pos & ts)
        print(f"  c{c}: 10k perms in {time.time()-t0:.0f}s")

        # Stats per test set
        for k, ts in test_sets.items():
            obs = obs_x[k]
            null = null_x[k]
            null_mean = float(null.mean()); null_sd = float(null.std())
            fold = obs / null_mean if null_mean > 0 else float("nan")
            z = (obs - null_mean) / (null_sd if null_sd > 0 else 1.0)
            p_one = float((null >= obs).sum() / n_perm_eff)
            rows.append({
                "component": c, "test_set": k,
                "n_top": top_n, "obs_overlap": int(obs),
                "matched_mean": null_mean, "matched_sd": null_sd,
                "fold_vs_matched": fold,
                "z_score": z,
                "empirical_p": p_one,
            })

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "permutation_all_components.csv", index=False)
    print(f"\nwrote {OUT/'permutation_all_components.csv'}")

    # Print pivoted z-score table
    print("\n=== z-score (matched permutation, one-sided) ===")
    z_piv = df.pivot(index="component", columns="test_set", values="z_score")
    z_piv = z_piv[["r_flip","C_flip","class_flip","HLA-B","HLA-C","MICB","MICA"]]
    print(z_piv.round(2).to_string())
    print("\n=== fold vs matched ===")
    fp = df.pivot(index="component", columns="test_set", values="fold_vs_matched")
    fp = fp[["r_flip","C_flip","class_flip","HLA-B","HLA-C","MICB","MICA"]]
    print(fp.round(2).to_string())
    print("\n=== empirical p (one-sided) ===")
    pp = df.pivot(index="component", columns="test_set", values="empirical_p")
    pp = pp[["r_flip","C_flip","class_flip","HLA-B","HLA-C","MICB","MICA"]]
    print(pp.round(4).to_string())

    # ===== Figure: 2-panel heatmap (z-score and fold vs matched) =====
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))

    # Panel A: z-score
    ax = axes[0]
    zmat = z_piv.values
    vmax_z = max(8, np.nanmax(np.abs(zmat)))
    im = ax.imshow(zmat, aspect="auto", cmap="RdBu_r",
                    vmin=-vmax_z, vmax=vmax_z, interpolation="nearest")
    ax.set_xticks(range(z_piv.shape[1])); ax.set_xticklabels(z_piv.columns, rotation=30, ha="right")
    ax.set_yticks(range(K)); ax.set_yticklabels([f"c{c}" for c in range(K)])
    for i in range(K):
        for j in range(z_piv.shape[1]):
            v = zmat[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.1f}", ha="center", va="center",
                         fontsize=7.5, color="#fff" if abs(v) > 4 else "#111")
    fig.colorbar(im, ax=ax, fraction=0.04, label="z-score (matched perm)")
    ax.set_title("A. z-score (matched permutation, 10,000 perms)", loc="left")

    # Panel B: fold vs matched
    ax = axes[1]
    fmat = fp.values
    im = ax.imshow(fmat, aspect="auto", cmap="RdBu_r",
                    vmin=0.5, vmax=1.5, interpolation="nearest")
    ax.set_xticks(range(fp.shape[1])); ax.set_xticklabels(fp.columns, rotation=30, ha="right")
    ax.set_yticks(range(K)); ax.set_yticklabels([f"c{c}" for c in range(K)])
    for i in range(K):
        for j in range(fp.shape[1]):
            v = fmat[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                         fontsize=7.5, color="#fff" if abs(v - 1) > 0.3 else "#111")
    fig.colorbar(im, ax=ax, fraction=0.04, label="fold vs matched mean")
    ax.set_title("B. fold vs matched mean", loc="left")

    fig.suptitle(
        "MAF×distance matched permutation (k=8, top 5%, 10,000 perms)\n"
        "Bins: MAF [0, 0.01, 0.05, 0.10, 0.20, 0.30, 1.0] × distance [0, 25, 50, 100, 150, 200, 260] kb",
        fontsize=10.5, y=1.04,
    )
    fig.tight_layout()
    out = OUT / "figure_permutation_all_components"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
