"""Supplementary Figure S7 — Reference-panel replication and NMF extension
of the COMT Val158Met (rs4680) flip-flop example (Lin et al. AJHG 2007).

Three panels:

  (A)  Signed r and asymmetry C per 1000G sub-population, for
       rs4680 (REF=G=Val) × rs2097603 (chr22:19,928,092 G>A, ALT=A).
  (B)  Component mixture (NMF k = 4) conditional on rs4680 allele —
       Val carrier vs Met carrier, per sub-population.
  (C)  Population-level COMT NMF mixture heatmap (26 pops × 4 components).

Inputs (existing intermediate tables from the upstream COMT analysis;
configure local path via the COMT_ANALYSIS_DIR environment variable,
default `data/comt_flipflop_nmf/`):
  $COMT_ANALYSIS_DIR/tables/comt_pair_metrics_by_pop.tsv
  $COMT_ANALYSIS_DIR/tables/comt_component_mixture_by_rs4680_status.tsv
  $COMT_ANALYSIS_DIR/tables/comt_component_mixture_by_pop.tsv

Outputs (mica_flipflop supplementary directory):
  supplementary/figures/figure_S4_comt_flipflop.pdf
  supplementary/figures/figure_S4_comt_flipflop.png
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

from figure_style import apply_style, cm_to_inch

REPO = Path(__file__).resolve().parent
COMT = Path(os.environ.get("COMT_ANALYSIS_DIR", "data/comt_flipflop_nmf"))

PAIR_TSV = COMT / "tables/comt_pair_metrics_by_pop.tsv"
RS4680_TSV = COMT / "tables/comt_component_mixture_by_rs4680_status.tsv"
POP_TSV = COMT / "tables/comt_component_mixture_by_pop.tsv"

OUT_DIR = REPO / "supplementary/figures"

SUPER_ORDER = ["AFR", "AMR", "EAS", "EUR", "SAS"]
SUPER_POP_MUTED = {
    "AFR": "#BFC4CC",
    "AMR": "#A6ACB5",
    "EAS": "#888E97",
    "EUR": "#5E6470",
    "SAS": "#3C414B",
}

COMP_COLORS = {
    "comp1": "#A6611A",
    "comp2": "#DDAA33",
    "comp3": "#6BAA3A",
    "comp4": "#2B7BBA",
}
COMP_LABELS = {
    "comp1": "comp1",
    "comp2": "comp2",
    "comp3": "comp3",
    "comp4": "comp4",
}


def _order_pops(df: pd.DataFrame) -> list[str]:
    pops = []
    for sp in SUPER_ORDER:
        sub = df[df["super_pop"] == sp]
        pops.extend(sorted(sub["pop"].unique()))
    return pops


def _draw_pop_bar(ax, x, vals, *, color_by_sign=False, facecolor=None,
                   ylabel: str = "", title: str = "", baseline: float = 0.0):
    width = 0.78
    if color_by_sign:
        colors = ["#2B7BBA" if v >= 0 else "#A6611A" for v in vals]
    else:
        colors = [facecolor] * len(vals)
    ax.bar(x, vals, width=width, color=colors,
            edgecolor="white", linewidth=0.3)
    ax.axhline(baseline, color="black", lw=0.6)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=9.5, fontweight="bold")
    ax.margins(x=0.005)


def _super_pop_strip(ax, pops, super_of, y=-0.15, h=0.045):
    for i, p in enumerate(pops):
        sp = super_of[p]
        ax.add_patch(plt.Rectangle((i - 0.4, y), 0.8, h,
                                     transform=ax.get_xaxis_transform(),
                                     clip_on=False,
                                     color=SUPER_POP_MUTED[sp], lw=0,
                                     alpha=0.9))


def panel_a(ax_r, ax_c, pair: pd.DataFrame, pops, super_of):
    pair = pair.set_index("pop").loc[pops]
    x = np.arange(len(pops))

    _draw_pop_bar(
        ax_r, x, pair["r"].values, color_by_sign=True,
        ylabel="signed r\n(rs4680 × rs2097603)",
        title="(A)  Population-conditional signed LD and asymmetry C")
    ax_r.set_ylim(-0.7, 0.7)
    ax_r.set_xticks(x)
    ax_r.set_xticklabels([])

    _draw_pop_bar(
        ax_c, x, pair["C"].values, color_by_sign=True,
        ylabel="asymmetry C\n(M + C decomposition)")
    cmax = max(0.05, np.abs(pair["C"]).max() * 1.15)
    ax_c.set_ylim(-cmax, cmax)
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(pops, rotation=90, fontsize=7)
    _super_pop_strip(ax_c, pops, super_of, y=-0.55, h=0.06)


def panel_b(ax, rs4680: pd.DataFrame, pops, super_of):
    val = rs4680[rs4680["status"] == "rs4680_Val_REF_carrier"].set_index("pop")
    met = rs4680[rs4680["status"] == "rs4680_Met_ALT_carrier"].set_index("pop")
    freq_cols = ["freq_comp1", "freq_comp2", "freq_comp3", "freq_comp4"]
    comps = ["comp1", "comp2", "comp3", "comp4"]

    val = val.loc[pops, freq_cols].copy()
    met = met.loc[pops, freq_cols].copy()
    val.columns = comps
    met.columns = comps
    val = val.div(val.sum(axis=1), axis=0)
    met = met.div(met.sum(axis=1), axis=0)

    n = len(pops)
    x = np.arange(n)
    bar_w = 0.36
    gap = 0.02
    x_val = x - (bar_w / 2 + gap / 2)
    x_met = x + (bar_w / 2 + gap / 2)

    for src, xpos, label in [(val, x_val, "Val"), (met, x_met, "Met")]:
        bottoms = np.zeros(n)
        for c in comps:
            v = src[c].values
            ax.bar(xpos, v, bottom=bottoms, width=bar_w,
                    color=COMP_COLORS[c], edgecolor="white", linewidth=0.3)
            bottoms = bottoms + v

    for xi in x_val:
        ax.text(xi, 1.02, "V", ha="center", va="bottom", fontsize=6,
                  color="#444444")
    for xi in x_met:
        ax.text(xi, 1.02, "M", ha="center", va="bottom", fontsize=6,
                  color="#444444")
    ax.set_xticks(x)
    ax.set_xticklabels(pops, rotation=90, fontsize=7)
    ax.set_ylim(0, 1.10)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("Component fraction\n(row-normalized)")
    ax.set_title(
        "(B)  k = 4 NMF component mixture conditional on rs4680 allele "
        "(V = Val/REF carrier, M = Met/ALT carrier)",
        loc="left", fontsize=9.5, fontweight="bold")
    _super_pop_strip(ax, pops, super_of, y=-0.30, h=0.045)


def panel_c(ax, pop_mix: pd.DataFrame, pops, super_of):
    comps = ["freq_comp1", "freq_comp2", "freq_comp3", "freq_comp4"]
    M = pop_mix.set_index("pop").loc[pops, comps].values
    im = ax.imshow(M.T, cmap="magma_r", aspect="auto",
                    vmin=0.0, vmax=max(0.60, M.max()))
    ax.set_xticks(np.arange(len(pops)))
    ax.set_xticklabels(pops, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(4))
    ax.set_yticklabels(["comp1", "comp2", "comp3", "comp4"], fontsize=8)
    ax.set_title(
        "(C)  Population-level COMT NMF mixture (k = 4, all 5,008 haplotypes)",
        loc="left", fontsize=9.5, fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, fraction=0.020, pad=0.012)
    cbar.set_label("mean component freq.", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)

    _super_pop_strip(ax, pops, super_of, y=-0.42, h=0.07)


def main():
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pair = pd.read_csv(PAIR_TSV, sep="\t")
    rs4680 = pd.read_csv(RS4680_TSV, sep="\t")
    pop_mix = pd.read_csv(POP_TSV, sep="\t")

    pops = _order_pops(pair)
    super_of = dict(zip(pair["pop"], pair["super_pop"]))

    fig = plt.figure(figsize=(cm_to_inch(22.0), cm_to_inch(22.0)))
    gs = GridSpec(4, 1, figure=fig,
                   height_ratios=[1.0, 1.0, 1.55, 1.0],
                   hspace=0.55)
    ax_a_r = fig.add_subplot(gs[0])
    ax_a_c = fig.add_subplot(gs[1], sharex=ax_a_r)
    ax_b = fig.add_subplot(gs[2])
    ax_c = fig.add_subplot(gs[3])

    panel_a(ax_a_r, ax_a_c, pair, pops, super_of)
    panel_b(ax_b, rs4680, pops, super_of)
    panel_c(ax_c, pop_mix, pops, super_of)

    comp_handles = [Patch(facecolor=COMP_COLORS[c], edgecolor="white",
                            label=COMP_LABELS[c]) for c in COMP_COLORS]
    sp_handles = [Patch(facecolor=SUPER_POP_MUTED[sp], label=sp)
                    for sp in SUPER_ORDER]
    sign_handles = [Patch(facecolor="#2B7BBA", label="r > 0 / C > 0"),
                      Patch(facecolor="#A6611A", label="r < 0 / C < 0")]

    leg1 = fig.legend(handles=sign_handles, loc="upper right",
                        bbox_to_anchor=(0.985, 0.985), frameon=False,
                        fontsize=7.5, title="Panel A sign",
                        title_fontsize=8, handlelength=1.2)
    leg1._legend_box.align = "left"
    leg2 = fig.legend(handles=comp_handles, loc="upper right",
                        bbox_to_anchor=(0.985, 0.62), frameon=False,
                        fontsize=7.5, title="NMF component (Panels B–C)",
                        title_fontsize=8, handlelength=1.2)
    leg2._legend_box.align = "left"
    leg3 = fig.legend(handles=sp_handles, loc="upper right",
                        bbox_to_anchor=(0.985, 0.36), frameon=False,
                        fontsize=7.5, title="Super-population",
                        title_fontsize=8, handlelength=1.2)
    leg3._legend_box.align = "left"

    fig.suptitle(
        "Supplementary Figure S7.  Reference-panel replication and NMF "
        "extension of the COMT Val158Met (rs4680) flip-flop example",
        y=0.985, va="top", wrap=True, fontsize=10.5, fontweight="bold")
    fig.tight_layout(rect=[0.015, 0.015, 0.83, 0.965])

    for ext in ("pdf", "png"):
        path = OUT_DIR / f"figure_S4_comt_flipflop.{ext}"
        fig.savefig(path, bbox_inches="tight",
                    dpi=300 if ext == "png" else None)
        print(f"wrote {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
