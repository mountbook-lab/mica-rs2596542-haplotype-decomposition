"""Threshold sensitivity (k=8 fixed, vary top % of NMF SNV signature).

Top thresholds: 1%, 2.5%, 5%, 10%, 15%
Per-component: r-flip, C-flip, class-flip, eQTL fold for MICA/HLA-B/HLA-C/MICB.

Output:
  results/threshold_sensitivity.csv
  results/figure_threshold_sensitivity.{pdf,png}
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
EQTL_CSV = REPO / "results/eqtl_targets_per_gene.csv"
REGION_CSV = REPO / "results/region_per_partner_26.csv"
OUT = REPO / "results"

THRESHOLDS = [0.01, 0.025, 0.05, 0.10, 0.15]
SUPER_POPS = ["AFR", "EUR", "EAS", "SAS", "AMR"]


def main():
    H = pd.read_parquet(H_PARQUET)
    snv_pos = np.array([int(c) for c in H.columns])
    h_mat = H.values
    K = h_mat.shape[0]
    n = len(snv_pos)
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

    rows = []
    for thr in THRESHOLDS:
        top_n = max(1, int(thr * n))
        for c in range(K):
            h = h_mat[c]
            top_idx = np.argsort(h)[-top_n:]
            top_set = set(int(p) for p in snv_pos[top_idx])
            row = {"threshold": thr, "component": c, "n_top": top_n}
            for label, ts in [("r_flip", rflip_set), ("C_flip", cflip_set),
                                ("class_flip", classflip_set),
                                ("MICA", eqtl_sets.get("MICA", set())),
                                ("MICB", eqtl_sets.get("MICB", set())),
                                ("HLA-B", eqtl_sets.get("HLA-B", set())),
                                ("HLA-C", eqtl_sets.get("HLA-C", set()))]:
                K_succ = len(ts); N = len(universe); n_d = top_n
                x = len(top_set & ts)
                if K_succ == 0:
                    fold = float("nan"); p = float("nan")
                else:
                    fold = x / (n_d * K_succ / N) if K_succ > 0 else float("nan")
                    p = float(hypergeom.sf(x - 1, N, K_succ, n_d))
                row[f"{label}_fold"] = fold; row[f"{label}_p"] = p
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "threshold_sensitivity.csv", index=False)

    print("=== Threshold sensitivity (fold enrichment) ===")
    cols = ["threshold","component","r_flip_fold","C_flip_fold","class_flip_fold",
             "MICA_fold","HLA-B_fold","HLA-C_fold","MICB_fold"]
    print(df[cols].to_string(index=False, formatters={
        "threshold":"{:.3f}".format,
        "r_flip_fold":"{:.2f}".format, "C_flip_fold":"{:.2f}".format,
        "class_flip_fold":"{:.2f}".format,
        "MICA_fold":"{:.2f}".format, "HLA-B_fold":"{:.2f}".format,
        "HLA-C_fold":"{:.2f}".format, "MICB_fold":"{:.2f}".format,
    }))

    # ===== Figure: line plot per component (focus on c5, c4, c6, c1) =====
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    fig, axes = plt.subplots(2, 4, figsize=(16, 6.5),
                              sharex=True, sharey=False)
    metrics = [("r_flip_fold", "r-flip fold"),
               ("C_flip_fold", "C-flip fold"),
               ("class_flip_fold", "class-flip fold"),
               ("MICA_fold", "MICA eQTL fold"),
               ("HLA-B_fold", "HLA-B eQTL fold"),
               ("HLA-C_fold", "HLA-C eQTL fold"),
               ("MICB_fold", "MICB eQTL fold"),
               (None, "(legend)")]
    palette = {0:"#888", 1:"#5a4fcf", 2:"#b5722d", 3:"#3070b8",
                4:"#3a8c5a", 5:"#d65454", 6:"#7e3da7", 7:"#1f4e8c"}
    style = {0:"-.", 1:"--", 2:"-", 3:":", 4:"-", 5:"-", 6:"-", 7:":"}
    for ax, (metric, title) in zip(axes.flat, metrics):
        if metric is None:
            ax.axis("off")
            for c in range(K):
                ax.plot([], [], color=palette[c], lw=1.5, ls=style[c],
                         label=f"c{c}", marker="o")
            ax.legend(loc="center", fontsize=10, ncol=2)
            continue
        for c in range(K):
            sub = df[df["component"] == c].sort_values("threshold")
            ax.plot(sub["threshold"]*100, sub[metric],
                     color=palette[c], lw=1.4, ls=style[c],
                     marker="o", markersize=4, label=f"c{c}",
                     alpha=0.85)
        ax.axhline(1.0, color="#222", lw=0.5, ls="--")
        ax.set_xlabel("top-% threshold")
        ax.set_ylabel(title)
        ax.set_title(title, loc="left")
        ax.grid(linestyle=":", alpha=0.3)
        ax.set_xticks([1, 2.5, 5, 10, 15])
    fig.suptitle(
        "Threshold sensitivity (NMF k=8, top % of SNV signature) — does the dichotomy hold across thresholds?",
        fontsize=11, y=1.01,
    )
    fig.tight_layout()
    out = OUT / "figure_threshold_sensitivity"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
