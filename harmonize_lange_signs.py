"""Add Lange-harmonized r_signed column to per_pop_topology CSVs and
re-render the Lange-replication figure with GRCh37 coordinates.

Convention table (locked, see rs2596542_allele_strand_build_mapping.md §4b):

  SNP          GRCh37 pos    dbSNP REF/ALT    Lange Fig 1A coded   sign_factor
  rs2596542    31,366,595    C / T            T  (anchor)           +1
  rs2244546    31,435,833    C / G            G                     +1
  rs9275572    32,678,999    A / G            A                     -1

For partner P, harmonized r = our r_signed × sign_factor[P], because
flipping the partner-coded allele negates the covariance.

Outputs:
  results/per_pop_topology_lange_coded.csv
  results/per_pop_topology_with_MC_lange_coded.csv
  results/figure_lange_replication_grch37.{pdf,png}
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
OUT = REPO / "results"

# Anchor and per-partner harmonization factors
PARTNER_INFO = {
    "rs2244546": {
        "grch37_pos": 31_435_833,
        "dbsnp_ref": "C",
        "dbsnp_alt": "G",
        "lange_coded": "G",
        "sign_factor": +1,
        "lange_jpt": +0.59,
        "lange_ceu": -0.19,
    },
    "rs9275572": {
        "grch37_pos": 32_678_999,
        "dbsnp_ref": "A",
        "dbsnp_alt": "G",
        "lange_coded": "A",
        "sign_factor": -1,
        "lange_jpt": +0.27,
        "lange_ceu": -0.12,
    },
}
ANCHOR_LABEL = "rs2596542 (chr6:31,366,595 GRCh37; coded T)"

CLASS_ORDER = ["DISJOINT", "STRICT_A_IN_B", "STRICT_B_IN_A",
                "PARTIAL", "IDENTICAL"]
CLASS_Y = {c: i for i, c in enumerate(CLASS_ORDER)}
SUPER_COLOR = {"EAS": "#5a4fcf", "EUR": "#3a8c5a"}


def add_harmonized_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add Lange-coding metadata + harmonized r_signed."""
    factor_map = {p: info["sign_factor"] for p, info in PARTNER_INFO.items()}
    pos_map    = {p: info["grch37_pos"]   for p, info in PARTNER_INFO.items()}
    coded_map  = {p: info["lange_coded"]  for p, info in PARTNER_INFO.items()}

    df = df.copy()
    df["partner_pos_grch37"] = df["partner"].map(pos_map)
    df["partner_coded_lange"] = df["partner"].map(coded_map)
    df["partner_sign_factor"] = df["partner"].map(factor_map)
    df["r_signed_lange_coded"] = (df["r_signed"].astype(float)
                                    * df["partner_sign_factor"])
    return df


def main():
    sys.stdout.reconfigure(line_buffering=True)
    OUT.mkdir(parents=True, exist_ok=True)

    # =============== (a) Augment CSVs ===============
    for src_name in ["per_pop_topology.csv", "per_pop_topology_with_MC.csv"]:
        src = OUT / src_name
        if not src.exists():
            print(f"  skip {src_name} (not found)"); continue
        df = pd.read_csv(src)
        df = add_harmonized_columns(df)
        out = OUT / src_name.replace(".csv", "_lange_coded.csv")
        df.to_csv(out, index=False)
        print(f"wrote {out}")
        # Quick sanity: print harmonized values for JPT, CEU
        sub = df[df["pop"].isin(["JPT", "CEU"])]
        print(sub[["partner", "pop", "super_pop", "r_signed",
                    "partner_sign_factor", "r_signed_lange_coded"]]
                  .to_string(index=False))
        print()

    # =============== (b) Re-render Lange-replication figure ===============
    df = pd.read_csv(OUT / "per_pop_topology_lange_coded.csv")

    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.0), sharey=True)

    rng = np.random.default_rng(0)
    for ax, (pname, info) in zip(axes, PARTNER_INFO.items()):
        sub = df[df["partner"] == pname].copy()
        sub["y"] = sub["topology_class"].map(CLASS_Y).astype(float)
        sub["y_jit"] = sub["y"] + rng.uniform(-0.18, 0.18, size=len(sub))

        # axis guides
        ax.axvline(0, color="#888", lw=0.7, ls="--")
        ax.axvline(+0.5, color="#b04545", lw=0.5, ls=":", alpha=0.5)
        ax.axvline(-0.5, color="#b04545", lw=0.5, ls=":", alpha=0.5)

        # Per-pop scatter (harmonized r on x-axis)
        for sup in ["EAS", "EUR"]:
            m = (sub["super_pop"] == sup) & (~sub["pop"].str.endswith("_pooled"))
            ax.scatter(sub.loc[m, "r_signed_lange_coded"], sub.loc[m, "y_jit"],
                        s=66, color=SUPER_COLOR[sup], edgecolor="#222",
                        linewidth=0.4, alpha=0.85, label=sup)
            mp = (sub["super_pop"] == sup) & sub["pop"].str.endswith("_pooled")
            ax.scatter(sub.loc[mp, "r_signed_lange_coded"], sub.loc[mp, "y"],
                        s=160, color=SUPER_COLOR[sup], edgecolor="#222",
                        linewidth=0.8, marker="D", alpha=0.95)
        # population labels
        for _, r in sub.iterrows():
            if r["pop"].endswith("_pooled"):
                continue
            ax.text(r["r_signed_lange_coded"] + 0.015, r["y_jit"], r["pop"],
                     fontsize=6.5, color="#333", va="center")

        # Mark Lange's published values as black ticks
        for pop_lab, sgn in [("JPT", info["lange_jpt"]),
                              ("CEU", info["lange_ceu"])]:
            ax.plot([sgn, sgn], [-0.65, -0.30], color="#222", lw=1.5)
            ax.text(sgn, -0.72, f"Lange\n{pop_lab}",
                     ha="center", va="top", fontsize=6.5, color="#222")

        ax.set_xlim(-0.7, 0.7)
        ax.set_yticks(range(len(CLASS_ORDER)))
        ax.set_yticklabels(CLASS_ORDER, fontsize=8)
        ax.set_xlabel(
            r"signed Pearson $r$  (Lange-coded: anchor=rs2596542-T, "
            f"partner={pname}-{info['lange_coded']})",
            fontsize=8.5,
        )
        sign_note = ("native dbSNP ALT coding = Lange coding"
                     if info["sign_factor"] == 1
                     else "harmonized × −1 (dbSNP ALT = G; Lange = A)")
        ax.set_title(
            f"rs2596542 × {pname}\n"
            f"chr6:{info['grch37_pos']:,} (GRCh37)  •  "
            f"Lange Fig 1A: r={info['lange_jpt']:+.2f} JPT, "
            f"{info['lange_ceu']:+.2f} CEU\n"
            f"({sign_note})",
            fontsize=9, loc="left",
        )
        ax.set_ylim(-1.0, len(CLASS_ORDER) - 0.5)
        ax.grid(axis="y", linestyle=":", alpha=0.4)
        if ax is axes[0]:
            ax.set_ylabel("topology class", fontsize=9)
        ax.legend(loc="upper right", fontsize=7.5, frameon=False)

    fig.suptitle(
        f"Lange 2013 Fig. 1A replication in 1000G phase 3, harmonized to Lange-coded alleles\n"
        f"Anchor: {ANCHOR_LABEL}.  Partner GRCh37 coordinates and coded alleles in panel titles.",
        fontsize=10.5, x=0.5, y=1.04,
    )
    fig.tight_layout()
    out = OUT / "figure_lange_replication_grch37"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"wrote {out}.pdf / .png")
    plt.close(fig)


if __name__ == "__main__":
    main()
