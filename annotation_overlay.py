"""Annotation overlay on NMF components (26-pop, k=8).

For each NMF component's top-5% signature SNVs, test enrichment for:
  - Gene-region membership (HCG27 / HLA-B / MICA / HCP5 / MICB / MICA-AS1)
  - Known HLA tag SNPs (paper 1 anchor set + Lange 2013 + Kumar 2011)
  - C-flip status (already done in decompose_branches_26.py, recomputed here
    for unified table)

Outputs:
  results/component_annotation_table.csv
  results/figure_component_annotation.{pdf,png}
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

REPO = Path(__file__).resolve().parent
H_PARQUET = REPO / "results/nmf_H_26_k8.parquet"
W_PARQUET = REPO / "results/nmf_W_26_k8.parquet"
REGION_CSV = REPO / "results/region_per_partner_26.csv"
OUT = REPO / "results"

# Gene / feature regions in the window (GRCh37)
GENE_REGIONS = [
    ("HCG27",     31_265_489, 31_272_009),
    ("HLA-B",     31_321_649, 31_324_989),
    ("MICA",      31_367_561, 31_383_090),
    ("MICA-AS1",  31_372_000, 31_381_000),
    ("HCP5",      31_430_950, 31_445_000),  # extended to lncRNA
    ("MICB",      31_462_658, 31_478_901),
    ("MCCD1",     31_496_488, 31_498_018),
]

# Known tag / anchor SNPs in window (paper 1 + Lange 2013)
TAG_SNPS = {
    "rs10484555 (HLA-B*15:02 tag)":     31_281_465,
    "rs2596542 (anchor, MICA up.)":     31_366_595,
    "rs2395029 (HLA-B*57:01 tag, HCP5)": 31_431_780,
    "rs2244546 (Lange 2013 HCP5 tag)":  31_435_833,
    "rs11509487 (MICB, paper 1)":       31_475_084,
}

NMF_TOP_FRAC = 0.05


def assign_gene(pos: int):
    """Return list of gene labels overlapping pos."""
    return [name for name, s, e in GENE_REGIONS if s <= pos <= e]


def main():
    H = pd.read_parquet(H_PARQUET)
    W = pd.read_parquet(W_PARQUET)
    region = pd.read_csv(REGION_CSV)
    cflip = set(region.loc[region["any_C_flip"] == True, "partner_pos"].astype(int))

    snv_pos = np.array([int(c) for c in H.columns])
    K = H.shape[0]
    h_mat = H.values
    print(f"NMF k={K}, n SNVs = {len(snv_pos)}, "
           f"n carriers = {len(W)}")

    # Gene overlap per SNV
    snv_genes = pd.DataFrame({
        "pos": snv_pos,
        "genes": ["+".join(assign_gene(p)) or "intergenic" for p in snv_pos],
    })

    # Universe = SNVs in NMF (already the partner set)
    universe = set(int(p) for p in snv_pos)
    region_uni = region[region["partner_pos"].isin(universe)]
    cflip_uni = set(region_uni.loc[region_uni["any_C_flip"] == True,
                                     "partner_pos"].astype(int))

    # Per-component top SNV set
    rows = []
    for c in range(K):
        h = h_mat[c]
        top_n = max(1, int(NMF_TOP_FRAC * len(h)))
        top_idx = np.argsort(h)[-top_n:]
        top_pos = snv_pos[top_idx]
        top_set = set(int(p) for p in top_pos)

        # Gene-region enrichment
        for gname, gs, ge in GENE_REGIONS:
            in_region_universe = (snv_pos >= gs) & (snv_pos <= ge)
            in_region_top = sum(gs <= p <= ge for p in top_pos)
            N = len(snv_pos); K_succ = int(in_region_universe.sum())
            n_draws = top_n; x = in_region_top
            if K_succ == 0:
                p_val = float("nan"); fold = float("nan")
            else:
                p_val = float(hypergeom.sf(x - 1, N, K_succ, n_draws))
                expected = n_draws * K_succ / N
                fold = x / expected if expected > 0 else float("nan")
            rows.append({
                "component": c, "annotation_type": "gene_region",
                "annotation": gname,
                "n_universe": K_succ, "n_top": x,
                "expected": n_draws * K_succ / N,
                "fold_enrichment": fold, "hyper_p": p_val,
            })

        # C-flip enrichment
        x_cf = len(top_set & cflip_uni)
        K_cf = len(cflip_uni); N = len(universe); n_draws = top_n
        p_cf = float(hypergeom.sf(x_cf - 1, N, K_cf, n_draws))
        rows.append({
            "component": c, "annotation_type": "C_flip",
            "annotation": "any_C_flip",
            "n_universe": K_cf, "n_top": x_cf,
            "expected": n_draws * K_cf / N,
            "fold_enrichment": x_cf / (n_draws * K_cf / N) if K_cf > 0 else float("nan"),
            "hyper_p": p_cf,
        })

        # Known tag SNPs presence (binary check, no enrichment test)
        for tag_name, tag_pos in TAG_SNPS.items():
            in_top = tag_pos in top_set
            in_universe = tag_pos in universe
            rows.append({
                "component": c, "annotation_type": "tag_snp",
                "annotation": tag_name,
                "n_universe": int(in_universe), "n_top": int(in_top),
                "expected": float("nan"), "fold_enrichment": float("nan"),
                "hyper_p": float("nan"),
            })
    annot = pd.DataFrame(rows)
    annot.to_csv(OUT / "component_annotation_table.csv", index=False)

    # Print pivoted summary
    print("\n=== Gene-region fold enrichment (per component, top 5% SNVs) ===")
    g = annot[annot["annotation_type"] == "gene_region"].pivot(
        index="component", columns="annotation",
        values="fold_enrichment").round(2)
    print(g.to_string())

    print("\n=== Gene-region hypergeom p-value (-log10) ===")
    p = annot[annot["annotation_type"] == "gene_region"].pivot(
        index="component", columns="annotation", values="hyper_p")
    print((-np.log10(p.replace(0, 1e-300))).round(1).to_string())

    print("\n=== Known tag SNPs in top 5% of each component (1=present) ===")
    t = annot[annot["annotation_type"] == "tag_snp"].pivot(
        index="component", columns="annotation", values="n_top")
    print(t.to_string())

    print("\n=== C-flip fold enrichment per component ===")
    cf = annot[annot["annotation_type"] == "C_flip"][["component","fold_enrichment","hyper_p"]]
    print(cf.to_string(index=False, formatters={
        "fold_enrichment": "{:.2f}".format, "hyper_p": "{:.2e}".format,
    }))

    # Figure: heatmap of gene-region fold enrichment + C-flip + tag presence
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
    })
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4),
                              gridspec_kw={"width_ratios": [3, 1.5, 2.4]})

    # Panel A: gene-region fold enrichment heatmap
    ax = axes[0]
    g_vals = g.values
    im = ax.imshow(g_vals, aspect="auto", cmap="RdBu_r",
                    vmin=0, vmax=4, interpolation="nearest")
    ax.set_xticks(range(g.shape[1])); ax.set_xticklabels(g.columns,
                                                          rotation=30, ha="right")
    ax.set_yticks(range(g.shape[0])); ax.set_yticklabels([f"c{c}" for c in g.index])
    for i in range(g.shape[0]):
        for j in range(g.shape[1]):
            v = g_vals[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                         fontsize=7, color="#111" if v < 2.5 else "#fff")
    ax.set_title("A. Gene-region fold enrichment (top 5% SNVs)", loc="left")
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("fold", fontsize=8)

    # Panel B: C-flip fold enrichment
    ax = axes[1]
    cf_vals = cf.set_index("component")["fold_enrichment"].values
    bars = ax.barh(range(len(cf_vals)), cf_vals,
                    color=["#d65454" if v >= 1.1 else "#888"
                           for v in cf_vals],
                    edgecolor="#222", linewidth=0.5)
    ax.axvline(1.0, color="#222", lw=0.5, ls="--")
    ax.set_yticks(range(len(cf_vals))); ax.set_yticklabels([f"c{c}" for c in cf.set_index("component").index])
    ax.invert_yaxis()
    ax.set_xlabel("fold enrichment (any-pair C-flip)")
    ax.set_title("B. C-flip enrichment", loc="left")
    for i, v in enumerate(cf_vals):
        ax.text(v + 0.02, i, f"{v:.2f}", va="center", fontsize=7)

    # Panel C: tag SNP presence matrix
    ax = axes[2]
    t_vals = t.values.astype(float)
    im = ax.imshow(t_vals, aspect="auto", cmap="Greys",
                    vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(t.shape[1]))
    ax.set_xticklabels([c.split(" ")[0] for c in t.columns],
                        rotation=30, ha="right", fontsize=7.5)
    ax.set_yticks(range(t.shape[0]))
    ax.set_yticklabels([f"c{c}" for c in t.index])
    for i in range(t.shape[0]):
        for j in range(t.shape[1]):
            v = t_vals[i, j]
            if v == 1:
                ax.text(j, i, "✓", ha="center", va="center",
                         fontsize=10, color="#1f4e8c")
    ax.set_title("C. Known tag SNPs in component top 5%", loc="left")

    fig.suptitle(
        f"NMF component annotation overlay (k=8, 26-pop, top {NMF_TOP_FRAC*100:.0f}% SNVs)",
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    out = OUT / "figure_component_annotation"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
