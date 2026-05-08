"""Render exploration figure: signed r vs topology class.

Two-panel:
  Left  -- rs2596542 x rs2244546 (Lange Fig 1A example)
  Right -- rs2596542 x rs9275572 (Lange Fig 1A example)

Each panel: x = signed Pearson r, y = topology class category.
Colour by super-pop (EAS purple / EUR green); marker shape by pop.

Output:
  figure_signed_r_vs_topology.{pdf,png}
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
CSV = REPO / "results/per_pop_topology.csv"

CLASS_ORDER = ["DISJOINT", "STRICT_A_IN_B", "STRICT_B_IN_A",
                "PARTIAL", "IDENTICAL"]
CLASS_Y = {c: i for i, c in enumerate(CLASS_ORDER)}

SUPER_COLOR = {"EAS": "#5a4fcf", "EUR": "#3a8c5a"}

LANGE = {
    "rs2244546": {"JPT": +0.59, "CEU": -0.19},
    "rs9275572": {"JPT": +0.27, "CEU": -0.12},
}


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
    })
    df = pd.read_csv(CSV)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)

    for ax, (pname, lange_signs) in zip(axes, LANGE.items()):
        sub = df[df["partner"] == pname].copy()
        sub["y"] = sub["topology_class"].map(CLASS_Y)

        # jitter y slightly per row to avoid overplotting
        rng = np.random.default_rng(0)
        sub["y_jit"] = sub["y"] + rng.uniform(-0.18, 0.18, size=len(sub))

        # vertical r=0 line
        ax.axvline(0, color="#888", lw=0.7, ls="--")
        # threshold lines
        ax.axvline(+0.5, color="#b04545", lw=0.6, ls=":", alpha=0.6)
        ax.axvline(-0.5, color="#b04545", lw=0.6, ls=":", alpha=0.6)

        # Plot per super-pop
        for sup in ["EAS", "EUR"]:
            m = (sub["super_pop"] == sup) & (~sub["pop"].str.endswith("_pooled"))
            ax.scatter(sub.loc[m, "r_signed"], sub.loc[m, "y_jit"],
                        s=60, color=SUPER_COLOR[sup], edgecolor="#222",
                        linewidth=0.4, alpha=0.85, label=sup)
            # Pooled markers (larger, no jitter)
            mp = (sub["super_pop"] == sup) & sub["pop"].str.endswith("_pooled")
            ax.scatter(sub.loc[mp, "r_signed"], sub.loc[mp, "y"],
                        s=140, color=SUPER_COLOR[sup], edgecolor="#222",
                        linewidth=0.8, marker="D", alpha=0.95)

        # Annotate population labels next to points
        for _, r in sub.iterrows():
            if r["pop"].endswith("_pooled"):
                continue
            ax.text(r["r_signed"] + 0.015, r["y_jit"], r["pop"],
                     fontsize=6.5, color="#333", va="center")

        # Mark Lange's reported values as vertical short bars
        for pop_lab, sgn in lange_signs.items():
            ax.plot([sgn, sgn], [-0.55, -0.30],
                     color="#222", lw=1.5)
            ax.text(sgn, -0.62, f"Lange\n{pop_lab}",
                     ha="center", va="top", fontsize=6.5, color="#222")

        ax.set_xlim(-0.7, 0.7)
        ax.set_yticks(range(len(CLASS_ORDER)))
        ax.set_yticklabels(CLASS_ORDER, fontsize=8)
        ax.set_xlabel(r"signed Pearson $r$  (anchor=rs2596542, "
                       r"REF/ALT coding = dbSNP forward)",
                       fontsize=8)
        ax.set_title(
            f"rs2596542 × {pname}\n"
            f"(Lange 2013 Fig 1A: "
            f"r={lange_signs['JPT']:+.2f} JPT, {lange_signs['CEU']:+.2f} CEU)",
            fontsize=9, loc="left",
        )
        ax.set_ylim(-0.9, len(CLASS_ORDER) - 0.5)
        ax.grid(axis="y", linestyle=":", alpha=0.4)
        if ax is axes[0]:
            ax.set_ylabel("topology class", fontsize=9)

        ax.legend(loc="upper right", fontsize=7.5)

    fig.suptitle(
        "Signed Pearson r vs carrier-set topology class — "
        "1000G phase 3, EAS + EUR sub-populations",
        fontsize=10, x=0.5, y=1.02,
    )
    fig.tight_layout()
    out = REPO / "results" / "figure_signed_r_vs_topology"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"wrote {out}.png and .pdf")


if __name__ == "__main__":
    main()
