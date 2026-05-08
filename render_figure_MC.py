"""Render M-C plane figure for Lange 2013 examples.

Two-panel: rs2596542 × rs2244546 and rs2596542 × rs9275572.
x = M (= p_A - p_B; frequency-difference component)
y = C (= Delta - M; structural / D-driven component)

Colour by super-pop. Sign(C) flip across populations is the cross-pop
flip captured by Delta=M+C decomposition.

Output: results/figure_MC_plane.{pdf,png}
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
CSV = REPO / "results/per_pop_topology_with_MC.csv"
SUPER_COLOR = {"EAS": "#5a4fcf", "EUR": "#3a8c5a"}


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
    })
    df = pd.read_csv(CSV)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)

    for ax, partner in zip(axes, ["rs2244546", "rs9275572"]):
        sub = df[df["partner"] == partner].copy()
        # Reference lines
        ax.axhline(0, color="#888", lw=0.7, ls="--")
        ax.axvline(0, color="#888", lw=0.7, ls="--")

        # Per-pop
        for sup in ["EAS", "EUR"]:
            m = (sub["super_pop"] == sup) & (~sub["pop"].str.endswith("_pooled"))
            ax.scatter(sub.loc[m, "M"], sub.loc[m, "C"],
                        s=70, color=SUPER_COLOR[sup], edgecolor="#222",
                        linewidth=0.4, alpha=0.9, label=sup)
            mp = (sub["super_pop"] == sup) & sub["pop"].str.endswith("_pooled")
            ax.scatter(sub.loc[mp, "M"], sub.loc[mp, "C"],
                        s=160, color=SUPER_COLOR[sup], edgecolor="#222",
                        linewidth=0.9, marker="D", alpha=0.95)

        # Population labels
        for _, r in sub.iterrows():
            if r["pop"].endswith("_pooled"): continue
            ax.text(r["M"] + 0.012, r["C"], r["pop"],
                     fontsize=6.5, color="#333", va="center")

        # Mark DISJOINT (topology class) with extra ring
        m_dj = (sub["topology_class"] == "DISJOINT")
        if m_dj.any():
            ax.scatter(sub.loc[m_dj, "M"], sub.loc[m_dj, "C"],
                        s=180, facecolors="none", edgecolor="#b04545",
                        linewidth=1.5, zorder=4)
            ax.text(sub.loc[m_dj, "M"].iloc[0],
                     sub.loc[m_dj, "C"].iloc[0] - 0.04,
                     "DISJOINT class", fontsize=7,
                     color="#b04545", ha="center", va="top")

        ax.set_xlabel(r"$M = p_A - p_B$  (frequency difference)", fontsize=8.5)
        if ax is axes[0]:
            ax.set_ylabel(r"$C = \Delta - M$  (structural / $D$-driven)",
                            fontsize=8.5)
        ax.set_title(f"rs2596542 × {partner}", fontsize=9.5, loc="left")
        ax.grid(linestyle=":", alpha=0.4)

        ax.legend(loc="upper right", fontsize=7.5)

    # Extras: text annotation explaining the flip pattern
    fig.text(
        0.5, -0.01,
        r"$\mathrm{sign}(C)$ flips perfectly between EAS and EUR for both pairs "
        r"(5/5 vs 5/5), even though topology class flips only for rs2244546 (CEU $\to$ DISJOINT). "
        r"$M$ sign is constant within each pair (frequency component does not flip).",
        ha="center", va="top", fontsize=8, color="#444",
    )
    fig.suptitle(
        r"$\Delta = M + C$ decomposition — Lange 2013 examples in 1000G phase 3",
        fontsize=10.5, y=1.00,
    )
    fig.tight_layout()
    out = REPO / "results" / "figure_MC_plane"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"wrote {out}.png and .pdf")


if __name__ == "__main__":
    main()
