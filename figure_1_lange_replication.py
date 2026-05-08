"""Figure 1 — Lange 2013 Fig. 1A replication in 1000G phase 3.

Two panels (rs2596542 × rs2244546 / rs9275572), signed r × 10 sub-pops.

  X-axis: 5 EAS (CHB, JPT, CHS, CDX, KHV) + 5 EUR (CEU, TSI, FIN, GBR, IBS)
  Y-axis: signed Pearson r in Lange coding
          (anchor=rs2596542-T; partner=G for rs2244546, A for rs9275572)
  Colored by super-pop; pooled super-pop shown as diamond marker on
  the right of each block.

  Lange Fig. 1A reported JPT / CEU values are overlaid as black ticks
  for direct visual comparison.

Input:
  results/per_pop_topology_lange_coded.csv  (built by harmonize_lange_signs.py)

Output:
  results/figure_1_lange_replication.{pdf,png}
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
SRC = REPO / "results/per_pop_topology_lange_coded.csv"
OUT = REPO / "results"

EAS_POPS = ["CHB", "JPT", "CHS", "CDX", "KHV"]
EUR_POPS = ["CEU", "TSI", "FIN", "GBR", "IBS"]
SUPER_COLOR = {"EAS": "#5a4fcf", "EUR": "#3a8c5a"}

PARTNER_INFO = {
    "rs2244546": {
        "grch37_pos": 31_435_833,
        "lange_coded": "G",
        "lange_jpt": +0.59,
        "lange_ceu": -0.19,
        "harmonization": "native dbSNP ALT = Lange coded G",
    },
    "rs9275572": {
        "grch37_pos": 32_678_999,
        "lange_coded": "A",
        "lange_jpt": +0.27,
        "lange_ceu": -0.12,
        "harmonization": r"$r$ × −1  (dbSNP ALT = G; Lange coded = A)",
    },
}


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
    })

    df = pd.read_csv(SRC)

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.8), sharex=True)

    pop_order = EAS_POPS + EUR_POPS
    pooled_x = len(pop_order) + 0.7   # diamond markers slightly to the right
    x_eas_pool = pooled_x
    x_eur_pool = pooled_x + 0.9

    for ax, (pname, info) in zip(axes, PARTNER_INFO.items()):
        sub = df[df["partner"] == pname].copy()
        sub = sub.set_index("pop")

        # Per-pop bars (5 EAS + 5 EUR)
        xs = np.arange(len(pop_order))
        rs = []
        cols = []
        for p in pop_order:
            if p in sub.index:
                r = float(sub.loc[p, "r_signed_lange_coded"])
                sup = sub.loc[p, "super_pop"]
                rs.append(r); cols.append(SUPER_COLOR[sup])
            else:
                rs.append(np.nan); cols.append("#888")
        rs = np.array(rs)
        bars = ax.bar(xs, rs, width=0.7, color=cols, edgecolor="#222",
                       linewidth=0.5, alpha=0.85)

        # Bar value labels
        for x, r in zip(xs, rs):
            if np.isnan(r): continue
            ax.text(x, r + (0.025 if r >= 0 else -0.025),
                     f"{r:+.2f}", ha="center",
                     va="bottom" if r >= 0 else "top",
                     fontsize=7.0, color="#222")

        # Pooled super-pop diamond markers
        for sup, xpool in [("EAS", x_eas_pool), ("EUR", x_eur_pool)]:
            pool_key = f"{sup}_pooled"
            if pool_key in sub.index:
                rp = float(sub.loc[pool_key, "r_signed_lange_coded"])
                ax.scatter([xpool], [rp], s=140, color=SUPER_COLOR[sup],
                            edgecolor="#222", linewidth=0.9, marker="D",
                            zorder=5)
                ax.text(xpool, rp + (0.04 if rp >= 0 else -0.04),
                         f"{rp:+.2f}", ha="center",
                         va="bottom" if rp >= 0 else "top",
                         fontsize=7.5, color="#222", fontweight="bold")
                ax.text(xpool, -0.85, f"{sup}\npooled",
                         ha="center", va="bottom", fontsize=7.0,
                         color=SUPER_COLOR[sup], fontweight="bold")

        # r=0 line
        ax.axhline(0, color="#222", lw=0.7)
        # EAS/EUR divider
        ax.axvline(len(EAS_POPS) - 0.5, color="#888", lw=0.6, ls=":")

        # Lange's reported JPT / CEU values as horizontal reference ticks
        # JPT is at xs[1], CEU at xs[5]
        jpt_x = pop_order.index("JPT")
        ceu_x = pop_order.index("CEU")
        for x_ref, lange_val, lab in [(jpt_x, info["lange_jpt"], "Lange JPT"),
                                        (ceu_x, info["lange_ceu"], "Lange CEU")]:
            ax.plot([x_ref - 0.45, x_ref + 0.45], [lange_val, lange_val],
                     color="#000", lw=1.6, ls="-", solid_capstyle="butt",
                     zorder=6)
            ax.text(x_ref, lange_val + (0.03 if lange_val >= 0 else -0.03),
                     f"{lab}\n{lange_val:+.2f}",
                     ha="center", va="bottom" if lange_val >= 0 else "top",
                     fontsize=6.5, color="#000",
                     fontweight="bold")

        # Y-axis cosmetics
        ax.set_ylabel(
            r"signed Pearson $r$" "\n"
            f"(coded: rs2596542-T × {pname}-{info['lange_coded']})",
            fontsize=8.5,
        )
        ax.set_ylim(-0.95, 0.85)
        ax.grid(axis="y", linestyle=":", alpha=0.35)

        # Title with GRCh37 position + harmonization
        ax.set_title(
            f"rs2596542 × {pname}  "
            f"(chr6:{info['grch37_pos']:,} GRCh37)  •  "
            f"Lange Fig. 1A: {info['lange_jpt']:+.2f} JPT / "
            f"{info['lange_ceu']:+.2f} CEU  •  {info['harmonization']}",
            loc="left", fontsize=9.5, pad=2,
        )

        # super-pop legend (only on top panel)
        if ax is axes[0]:
            from matplotlib.patches import Patch
            handles = [Patch(facecolor=SUPER_COLOR["EAS"], edgecolor="#222",
                                  label="EAS sub-pops"),
                        Patch(facecolor=SUPER_COLOR["EUR"], edgecolor="#222",
                                  label="EUR sub-pops")]
            ax.legend(handles=handles, loc="upper right", fontsize=7.5,
                      frameon=False)

    # X-axis (bottom panel only)
    axes[1].set_xticks(list(range(len(pop_order))) + [x_eas_pool, x_eur_pool])
    axes[1].set_xticklabels(pop_order + ["", ""], fontsize=8.5)
    axes[1].set_xlabel("1000 Genomes Phase 3 sub-population",
                       fontsize=9, labelpad=4)
    axes[0].set_xticks(list(range(len(pop_order))) + [x_eas_pool, x_eur_pool])
    axes[0].set_xticklabels([""] * (len(pop_order) + 2))

    fig.suptitle(
        "Figure 1.  Population-conditional sign of LD between rs2596542 (anchor) "
        "and two HCV-HCC tag candidates — replicating Lange et al. 2013 (J Hepatol) Fig. 1A "
        "in 1000G phase 3.\n"
        "Anchor: rs2596542  (chr6:31,366,595 GRCh37, coded T).  "
        "Each bar = one of 10 EAS/EUR sub-populations; diamond = pooled super-pop. "
        "Black ticks = Lange's published JPT/CEU values.",
        fontsize=10, x=0.5, y=1.02,
    )
    fig.tight_layout()
    out = OUT / "figure_1_lange_replication"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"wrote {out}.pdf / .png")
    plt.close(fig)


if __name__ == "__main__":
    main()
