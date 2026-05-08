"""Region-wide figure: M-C plane and spatial scan around rs2596542.

Two panels:
  Left  -- median(C) EAS vs median(C) EUR scatter for all partners.
           Color-coded by C-flip status (NO flip / EAS+→EUR- / EAS-→EUR+).
           DISJOINT-class partners overlaid as red rings.
  Right -- chromosomal position vs median(C) EUR - median(C) EAS
           (delta-C across super-pops). Anchor position marked.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

REPO = Path(__file__).resolve().parent
PARQUET = REPO / "results/region_per_pop.parquet"
ANCHOR_POS = 31_366_595


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
    })
    df = pd.read_parquet(PARQUET)
    EAS = ["CHB", "JPT", "CHS", "CDX", "KHV"]
    EUR = ["CEU", "TSI", "FIN", "GBR", "IBS"]
    df["super"] = df["pop"].apply(
        lambda p: "EAS" if p in EAS else ("EUR" if p in EUR else "pooled"))

    sub = df[df["super"].isin(["EAS","EUR"])].copy()
    g = sub.groupby(["partner_pos", "super"]).agg(
        median_C=("C", "median"),
        median_M=("M", "median"),
        median_r=("r_signed", "median"),
        # mode topology class (most common in 5 pops)
        mode_class=("topology_class",
                     lambda s: s.value_counts().idxmax()),
        any_disjoint=("topology_class", lambda s: (s == "DISJOINT").any()),
    ).reset_index()

    wide = g.pivot(index="partner_pos", columns="super",
                    values=["median_C","median_M","median_r","mode_class","any_disjoint"]).reset_index()
    wide.columns = ["_".join(filter(None, [a, b])) for a, b in wide.columns]
    wide = wide.rename(columns={"partner_pos_": "partner_pos"})
    wide = wide.dropna(subset=["median_C_EAS","median_C_EUR"]).copy()

    # C-flip status
    def flip_status(eas, eur):
        if eas > 0 and eur < 0: return "flip_pos_neg"
        if eas < 0 and eur > 0: return "flip_neg_pos"
        if eas > 0 and eur > 0: return "no_flip_pos"
        if eas < 0 and eur < 0: return "no_flip_neg"
        return "zero"
    wide["c_status"] = [flip_status(e, u) for e, u in
                         zip(wide["median_C_EAS"], wide["median_C_EUR"])]
    wide["delta_C"] = wide["median_C_EUR"] - wide["median_C_EAS"]
    wide["mode_flip"] = (wide["mode_class_EAS"] != wide["mode_class_EUR"])
    wide["any_disjoint_flip"] = (
        (wide["any_disjoint_EAS"] != wide["any_disjoint_EUR"])
    )

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))

    # ---- Panel A: median(C) EAS vs EUR scatter ----
    ax = axes[0]
    palette = {
        "no_flip_pos": "#4078a3",
        "no_flip_neg": "#4078a3",
        "flip_pos_neg": "#d65454",
        "flip_neg_pos": "#d65454",
        "zero": "#aaaaaa",
    }
    for status in ["no_flip_pos","no_flip_neg","zero","flip_pos_neg","flip_neg_pos"]:
        m = wide["c_status"] == status
        if not m.any(): continue
        ax.scatter(
            wide.loc[m, "median_C_EAS"], wide.loc[m, "median_C_EUR"],
            s=8, color=palette[status], alpha=0.5, edgecolor="none",
            label="C-flip" if status.startswith("flip") and status == "flip_pos_neg"
                  else ("no flip" if status == "no_flip_pos"
                        else (None if status != "flip_neg_pos" else None))
        )

    # Overlay class-flip (mode change) markers
    m = wide["mode_flip"]
    ax.scatter(
        wide.loc[m, "median_C_EAS"], wide.loc[m, "median_C_EUR"],
        s=22, facecolors="none", edgecolor="#222", linewidth=0.6,
        label=f"mode-class flip ({m.sum()})",
    )
    ax.axhline(0, color="#666", lw=0.6, ls="--")
    ax.axvline(0, color="#666", lw=0.6, ls="--")
    ax.plot([-0.5, 0.5], [-0.5, 0.5], color="#666", lw=0.5, alpha=0.5)
    ax.set_xlim(-0.5, 0.5); ax.set_ylim(-0.5, 0.5)
    ax.set_xlabel("median $C$ across EAS sub-pops"); ax.set_ylabel("median $C$ across EUR sub-pops")
    ax.set_title("A. $C$-flip across populations (anchor=rs2596542)", loc="left", fontsize=10)
    ax.grid(linestyle=":", alpha=0.4)

    # Quadrant counts
    n = len(wide)
    n_flip_pn = (wide["c_status"]=="flip_pos_neg").sum()
    n_flip_np = (wide["c_status"]=="flip_neg_pos").sum()
    n_no = ((wide["c_status"].str.startswith("no_flip"))).sum()
    ax.text(0.97, 0.05,
             f"$n$ = {n} partners\n"
             f"C-flip (EAS+/EUR−): {n_flip_pn} ({n_flip_pn/n*100:.1f}%)\n"
             f"C-flip (EAS−/EUR+): {n_flip_np} ({n_flip_np/n*100:.1f}%)\n"
             f"no flip:           {n_no} ({n_no/n*100:.1f}%)\n"
             f"mode-class flip:   {wide['mode_flip'].sum()} ({wide['mode_flip'].sum()/n*100:.1f}%)",
             transform=ax.transAxes, fontsize=7.5, ha="right", va="bottom",
             family="monospace",
             bbox=dict(boxstyle="round,pad=0.4", fc="#ffffff", ec="#888", lw=0.5))
    ax.legend(loc="upper left", fontsize=7.5)

    # ---- Panel B: spatial scan (position vs delta_C) ----
    ax = axes[1]
    color_b = wide["c_status"].map({
        "no_flip_pos": "#4078a3", "no_flip_neg": "#4078a3",
        "flip_pos_neg": "#d65454", "flip_neg_pos": "#d65454",
        "zero": "#aaaaaa"})
    ax.scatter(wide["partner_pos"]/1e6, wide["delta_C"],
                s=6, c=color_b, alpha=0.5, edgecolor="none")
    # Highlight class-flip partners
    m = wide["mode_flip"]
    ax.scatter(wide.loc[m, "partner_pos"]/1e6, wide.loc[m, "delta_C"],
                s=22, facecolors="none", edgecolor="#222", linewidth=0.6,
                label=f"mode-class flip ({m.sum()})")
    # Lange examples
    for pos, name in [(31_366_595, "rs2596542"),
                       (31_435_833, "rs2244546"),
                       (32_678_999, "rs9275572 (out of window)")]:
        if pos < 31_116_595 or pos > 31_616_595: continue
        ax.axvline(pos/1e6, color="#222", lw=0.5, ls=":", alpha=0.6)
        ax.text(pos/1e6, ax.get_ylim()[1] * 0.95, name,
                 fontsize=6.5, rotation=90, va="top", ha="right", color="#444")
    ax.axhline(0, color="#666", lw=0.6, ls="--")
    ax.set_xlabel("chromosomal position (Mb, GRCh37)")
    ax.set_ylabel(r"median $C_\mathrm{EUR}$ − median $C_\mathrm{EAS}$")
    ax.set_title("B. Spatial scan: $\\Delta C$ between super-populations", loc="left", fontsize=10)
    ax.grid(linestyle=":", alpha=0.3)
    ax.legend(loc="upper right", fontsize=7.5)

    fig.suptitle("Region-wide $\\Delta = M + C$ scan around rs2596542 (±250 kb)",
                  fontsize=11, y=1.01)
    fig.tight_layout()
    out = REPO / "results" / "figure_region_MC"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"wrote {out}.png and .pdf")

    # Stats summary
    print()
    print("=== Region-wide stats ===")
    print(f"  total partners: {len(wide)}")
    print(f"  C-flip rate:        {(wide['c_status'].str.startswith('flip')).mean()*100:.1f}%")
    print(f"  mode-class flip:    {wide['mode_flip'].mean()*100:.1f}%")
    print(f"  any-DISJOINT flip:  {wide['any_disjoint_flip'].mean()*100:.1f}%")
    # Among C-flip partners, fraction that also have class flip
    cflip = wide["c_status"].str.startswith("flip")
    print()
    print(f"  partners with C-flip: {cflip.sum()}")
    print(f"    of which mode-class also flips: {(cflip & wide['mode_flip']).sum()} "
          f"({(cflip & wide['mode_flip']).sum()/cflip.sum()*100:.1f}%)")
    print(f"    of which mode-class STAYS:      {(cflip & ~wide['mode_flip']).sum()} "
          f"({(cflip & ~wide['mode_flip']).sum()/cflip.sum()*100:.1f}%)")


if __name__ == "__main__":
    main()
