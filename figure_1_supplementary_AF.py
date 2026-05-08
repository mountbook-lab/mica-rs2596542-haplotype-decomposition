"""Figure 1 supplementary — allele-frequency comparison across populations.

Motivates the magnitude differences between 1000G phase 3 (this study) and
HapMap3 (Lange 2013). Same-direction signs but different magnitudes are
expected because:
  r = D / sqrt(p_a (1-p_a) p_b (1-p_b))
so when p_a or p_b shifts even slightly between cohorts, |r| can change
substantially even when D stays similar.

Three panels (each = 10 1000G phase 3 sub-pops, EAS vs EUR):

  Panel A: anchor — rs2596542-T   (chr6:31,366,595 GRCh37)
           Lange's reported values: 0.329 JPT, 0.284 CEU (2013 Results)
  Panel B: partner — rs2244546-G  (chr6:31,435,833 GRCh37)
  Panel C: partner — rs9275572-A  (chr6:32,678,999 GRCh37)
           note: stored p_b is for dbSNP ALT=G → freq(A) = 1 − p_b

Output:
  results/figure_1_supplementary_AF.{pdf,png}
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
SRC = REPO / "results/per_pop_topology_with_MC_lange_coded.csv"
OUT = REPO / "results"

EAS_POPS = ["CHB", "JPT", "CHS", "CDX", "KHV"]
EUR_POPS = ["CEU", "TSI", "FIN", "GBR", "IBS"]
ALL_POPS = EAS_POPS + EUR_POPS
SUPER_COLOR = {"EAS": "#5a4fcf", "EUR": "#3a8c5a"}

# Lange 2013 reported HapMap3 MAFs (only anchor reported in abstract/results)
LANGE_REPORTED = {
    "rs2596542-T": {"JPT": 0.329, "CEU": 0.284},
}


def af_per_pop(df: pd.DataFrame, partner: str, lange_flip: bool):
    """Return ordered list of (Lange-coded) allele frequencies across ALL_POPS
    plus pooled EAS / EUR.
    """
    sub = df[df["partner"] == partner].set_index("pop")
    out = {}
    for p in ALL_POPS + ["EAS_pooled", "EUR_pooled"]:
        if p in sub.index:
            pb = float(sub.loc[p, "p_b"])
            out[p] = (1.0 - pb) if lange_flip else pb
        else:
            out[p] = float("nan")
    return out


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
    })

    df = pd.read_csv(SRC)

    # Anchor p_a is the same for any partner (same anchor) — pull from rs2244546 row
    sub = df[df["partner"] == "rs2244546"].set_index("pop")
    anchor_af = {p: (float(sub.loc[p, "p_a"]) if p in sub.index else float("nan"))
                  for p in ALL_POPS + ["EAS_pooled", "EUR_pooled"]}

    rs2244546_af = af_per_pop(df, "rs2244546", lange_flip=False)
    rs9275572_af = af_per_pop(df, "rs9275572", lange_flip=True)

    panels = [
        ("A", "rs2596542-T  (anchor)\nchr6:31,366,595 GRCh37",
         anchor_af, LANGE_REPORTED["rs2596542-T"], "p(T) — anchor coded"),
        ("B", "rs2244546-G  (partner)\nchr6:31,435,833 GRCh37",
         rs2244546_af, None, "p(G) — Lange coded"),
        ("C", "rs9275572-A  (partner)\nchr6:32,678,999 GRCh37",
         rs9275572_af, None, "p(A) — Lange coded (= 1 − p(G))"),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(10.5, 8.5), sharex=True)

    pooled_x_eas = len(ALL_POPS) + 0.7
    pooled_x_eur = pooled_x_eas + 0.9

    for ax, (tag, title, af, lange_vals, ylab) in zip(axes, panels):
        xs = np.arange(len(ALL_POPS))
        rs = np.array([af[p] for p in ALL_POPS])
        cols = [SUPER_COLOR["EAS"] if p in EAS_POPS else SUPER_COLOR["EUR"]
                  for p in ALL_POPS]
        ax.bar(xs, rs, width=0.7, color=cols, edgecolor="#222",
                 linewidth=0.5, alpha=0.85)
        for x, r in zip(xs, rs):
            if np.isnan(r): continue
            ax.text(x, r + 0.012, f"{r:.2f}", ha="center", va="bottom",
                     fontsize=7.0, color="#222")

        # Pooled diamonds
        for sup, xpool in [("EAS", pooled_x_eas), ("EUR", pooled_x_eur)]:
            v = af.get(f"{sup}_pooled", float("nan"))
            if not np.isnan(v):
                ax.scatter([xpool], [v], s=140, color=SUPER_COLOR[sup],
                            edgecolor="#222", linewidth=0.9, marker="D",
                            zorder=5)
                ax.text(xpool, v + 0.018, f"{v:.2f}", ha="center", va="bottom",
                         fontsize=7.5, color="#222", fontweight="bold")
                ax.text(xpool, -0.04, f"{sup}\npooled",
                         ha="center", va="top", fontsize=7.0,
                         color=SUPER_COLOR[sup], fontweight="bold")

        # Lange reference horizontal segments at JPT / CEU bars (anchor only)
        if lange_vals is not None:
            for pop_lab, val in lange_vals.items():
                idx = ALL_POPS.index(pop_lab)
                ax.plot([idx - 0.45, idx + 0.45], [val, val],
                         color="#000", lw=1.6, ls="-", zorder=6)
                ax.text(idx, val + 0.03,
                         f"Lange {pop_lab}\n{val:.2f}",
                         ha="center", va="bottom", fontsize=6.5,
                         color="#000", fontweight="bold")

        # EAS / EUR divider
        ax.axvline(len(EAS_POPS) - 0.5, color="#888", lw=0.6, ls=":")

        ax.set_ylim(-0.10, 0.95)
        ax.set_ylabel(ylab, fontsize=8.5)
        ax.grid(axis="y", linestyle=":", alpha=0.35)
        ax.set_title(f"{tag}.  {title}", loc="left", fontsize=9.5, pad=2)

    # X-axis (bottom panel only)
    axes[-1].set_xticks(list(range(len(ALL_POPS))) + [pooled_x_eas, pooled_x_eur])
    axes[-1].set_xticklabels(ALL_POPS + ["", ""], fontsize=8.5)
    axes[-1].set_xlabel("1000 Genomes Phase 3 sub-population",
                         fontsize=9, labelpad=4)

    fig.suptitle(
        "Figure 1 — supplementary.  Allele-frequency context for the Lange Fig. 1A "
        "magnitude differences.\n"
        r"Recall $|r| = D / \sqrt{p_a(1-p_a)\,p_b(1-p_b)}$ — even small AF shifts "
        "between HapMap3 and 1000G can change magnitudes while preserving sign.\n"
        "Coordinates GRCh37; alleles in Lange 2013 coding (rs2596542-T, "
        "rs2244546-G, rs9275572-A).",
        fontsize=10, x=0.5, y=1.02,
    )
    fig.tight_layout()
    out = OUT / "figure_1_supplementary_AF"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"wrote {out}.pdf / .png")
    plt.close(fig)


if __name__ == "__main__":
    main()
