"""r-flip vs C-flip vs class-flip — comparison enrichment per NMF component.

Three flip definitions (all multi-super-pop, "any-pair flip" criterion):
  r_flip      : sign(median_r) flips for >=1 super-pop pair (10 pairs)
  C_flip      : sign(median_C) flips for >=1 super-pop pair  (already in CSV)
  class_flip  : mode_class differs across super-pops         (already in CSV)

For each NMF component (k=8, 26-pop):
  hypergeometric test of top 5% SNV overlap with each flip set.
Plus Venn-style overlap of the three flip sets across the entire universe.

Outputs:
  results/component_x_flip_type.csv
  results/figure_flip_type_comparison.{pdf,png}
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

REPO = Path(__file__).resolve().parent
H_PARQUET = REPO / "results/nmf_H_26_k8.parquet"
REGION_CSV = REPO / "results/region_per_partner_26.csv"
OUT = REPO / "results"

NMF_TOP_FRAC = 0.05
SUPER_POPS = ["AFR", "EUR", "EAS", "SAS", "AMR"]


def main():
    H = pd.read_parquet(H_PARQUET)
    snv_pos = np.array([int(c) for c in H.columns])
    K = H.shape[0]
    print(f"NMF k={K}, n SNVs = {len(snv_pos)}")

    region = pd.read_csv(REGION_CSV)
    print(f"region partners: {len(region)}")

    # ===== r_flip: same definition as C_flip but on median_r =====
    rmat = region[[f"median_r_{s}" for s in SUPER_POPS]].astype(float).values
    n_r_flips = np.zeros(len(region), dtype=int)
    n_pairs = np.zeros(len(region), dtype=int)
    for i, j in combinations(range(len(SUPER_POPS)), 2):
        ri = rmat[:, i]; rj = rmat[:, j]
        valid = ~(np.isnan(ri) | np.isnan(rj))
        n_pairs += valid.astype(int)
        flip = valid & ((np.sign(ri) * np.sign(rj)) < 0)
        n_r_flips += flip.astype(int)
    region["n_super_r_flips"] = n_r_flips
    region["any_r_flip"] = (n_r_flips >= 1)
    region["high_r_flip"] = (n_r_flips >= 4)

    # ===== Universes =====
    universe = set(int(p) for p in snv_pos)
    cflip_set = set(region.loc[region["any_C_flip"] == True, "partner_pos"].astype(int)) & universe
    rflip_set = set(region.loc[region["any_r_flip"] == True, "partner_pos"].astype(int)) & universe
    classflip_set = set(region.loc[region["any_class_flip"] == True, "partner_pos"].astype(int)) & universe

    print(f"\nUniverse size: {len(universe)}")
    print(f"  any C-flip:     {len(cflip_set)} ({len(cflip_set)/len(universe)*100:.1f}%)")
    print(f"  any r-flip:     {len(rflip_set)} ({len(rflip_set)/len(universe)*100:.1f}%)")
    print(f"  any class-flip: {len(classflip_set)} ({len(classflip_set)/len(universe)*100:.1f}%)")
    print()
    # Pairwise overlap (Jaccard)
    def jacc(a, b):
        if not a or not b: return 0.0
        return len(a & b) / len(a | b)
    print("Pairwise Jaccard:")
    print(f"  C ∩ r:     {len(cflip_set & rflip_set):5d}  Jaccard={jacc(cflip_set, rflip_set):.3f}")
    print(f"  C ∩ class: {len(cflip_set & classflip_set):5d}  Jaccard={jacc(cflip_set, classflip_set):.3f}")
    print(f"  r ∩ class: {len(rflip_set & classflip_set):5d}  Jaccard={jacc(rflip_set, classflip_set):.3f}")
    # Triple intersection
    triple = cflip_set & rflip_set & classflip_set
    print(f"  C ∩ r ∩ class: {len(triple)}")
    # C-only, r-only, etc.
    print(f"  C-only (not r):     {len(cflip_set - rflip_set)}")
    print(f"  r-only (not C):     {len(rflip_set - cflip_set)}")
    print(f"  C ∩ r:              {len(cflip_set & rflip_set)}")

    # ===== Per-component enrichment =====
    h_mat = H.values
    rows = []
    for c in range(K):
        h = h_mat[c]
        top_n = max(1, int(NMF_TOP_FRAC * len(h)))
        top_idx = np.argsort(h)[-top_n:]
        top_set = set(int(p) for p in snv_pos[top_idx])
        n_draws = len(top_set)
        N = len(universe)
        for label, fset in [("C_flip", cflip_set), ("r_flip", rflip_set),
                              ("class_flip", classflip_set)]:
            K_succ = len(fset)
            x = len(top_set & fset)
            if K_succ == 0:
                rows.append({"component": c, "flip_type": label,
                              "x": x, "expected": 0, "fold": float("nan"),
                              "p": float("nan")})
                continue
            p = float(hypergeom.sf(x - 1, N, K_succ, n_draws))
            exp = n_draws * K_succ / N
            rows.append({
                "component": c, "flip_type": label,
                "x": x, "expected": exp,
                "fold": x / exp if exp > 0 else float("nan"),
                "p": p,
            })
    enr = pd.DataFrame(rows)
    enr["neglog10_p"] = -np.log10(enr["p"].replace(0, 1e-300))
    enr.to_csv(OUT / "component_x_flip_type.csv", index=False)

    # Pivot for printing
    fold_pivot = enr.pivot(index="component", columns="flip_type", values="fold")
    p_pivot = enr.pivot(index="component", columns="flip_type", values="neglog10_p")
    print("\n=== Component × flip-type fold enrichment ===")
    print(fold_pivot[["r_flip","C_flip","class_flip"]].round(2).to_string())
    print("\n=== -log10 p ===")
    print(p_pivot[["r_flip","C_flip","class_flip"]].round(1).to_string())

    # ===== Figure: 4-panel =====
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
    })
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 4.6),
                              gridspec_kw={"width_ratios": [2, 2, 2, 2]})

    # Panels A-C: per-flip-type fold heatmap
    flip_types = [("r_flip", "A. r-flip"),
                   ("C_flip", "B. C-flip"),
                   ("class_flip", "C. class-flip")]
    fmat = fold_pivot[["r_flip","C_flip","class_flip"]].values
    pmat = p_pivot[["r_flip","C_flip","class_flip"]].values
    vmax = 2.0
    for ax_i, (col, title) in enumerate(flip_types):
        ax = axes[ax_i]
        vals = fold_pivot[col].values
        sig = -np.log10(enr[enr["flip_type"]==col].sort_values("component")["p"].values + 1e-300)
        bars = ax.barh(range(K), vals,
                        color=["#d65454" if (v >= 1.1 and s >= 1.3)
                               else "#888" for v, s in zip(vals, sig)],
                        edgecolor="#222", linewidth=0.5)
        ax.axvline(1.0, color="#222", lw=0.5, ls="--")
        ax.set_yticks(range(K)); ax.set_yticklabels([f"c{c}" for c in range(K)])
        ax.invert_yaxis()
        ax.set_xlabel(f"fold ({col})")
        ax.set_title(title, loc="left")
        ax.set_xlim(0, max(2.0, vals.max() * 1.05))
        for i, (v, s) in enumerate(zip(vals, sig)):
            ax.text(v + 0.02, i,
                     f"{v:.2f}\n(-log10p={s:.0f})",
                     va="center", fontsize=6.5)

    # Panel D: r-flip vs C-flip overlap of universe sets (Venn-style bars)
    ax = axes[3]
    universe_n = len(universe)
    only_C = len(cflip_set - rflip_set)
    only_r = len(rflip_set - cflip_set)
    both_Cr = len(cflip_set & rflip_set)
    only_class = len(classflip_set - cflip_set - rflip_set)
    triple = len(cflip_set & rflip_set & classflip_set)
    none = universe_n - len(cflip_set | rflip_set | classflip_set)
    cats = ["C only", "r only", "C ∩ r", "class\nonly", "C∩r∩class", "no flip"]
    counts = [only_C, only_r, both_Cr, only_class, triple, none]
    cols = ["#3a8c5a", "#5a4fcf", "#b04545", "#b5722d", "#7e3da7", "#aaaaaa"]
    ax.bar(cats, counts, color=cols, edgecolor="#222", linewidth=0.4)
    for i, c in enumerate(counts):
        ax.text(i, c + universe_n*0.005, f"{c}",
                 ha="center", va="bottom", fontsize=7)
    ax.set_ylabel("partner count")
    ax.set_title("D. Universe overlap of flip definitions", loc="left")
    ax.tick_params(axis="x", labelsize=7)

    fig.suptitle(
        "r-flip vs C-flip vs class-flip — enrichment per NMF component (k=8, 26-pop)",
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    out = OUT / "figure_flip_type_comparison"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
