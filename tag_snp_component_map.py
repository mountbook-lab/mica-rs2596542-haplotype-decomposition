"""Map known tag SNPs to NMF components.

Tag SNPs:
  rs2596542  31,366,595  (anchor — excluded from NMF, but quote its
                          neighborhood)
  rs2395029  31,431,780  (HCP5, HLA-B*57:01 tag — paper 1 case study)
  rs2244546  31,435,833  (HCP5 region — Lange 2013 better tag proposal)
  rs11509487 31,475,084  (MICB indel — paper 1 anchor; INDEL, may be filtered)

For each tag SNP that is in the NMF SNV set:
  - direct H[c, pos] for each component c
  - rank within component (where it sits relative to top-5% threshold)
  - nearest top-5% SNV (distance) per component

Output:
  results/tag_snp_component_map.csv
  results/figure_tag_snp_component.{pdf,png}
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
H_PARQUET = REPO / "results/nmf_H_26_k8.parquet"
OUT = REPO / "results"

TAG_SNPS = {
    "rs2596542 (anchor)":             31_366_595,
    "rs2395029 (HCP5, B*57:01 tag)":  31_431_780,
    "rs2244546 (Lange 2013 tag)":     31_435_833,
    "rs11509487 (MICB, paper 1)":     31_475_084,
}


def main():
    H = pd.read_parquet(H_PARQUET)
    snv_pos = np.array([int(c) for c in H.columns])
    h_mat = H.values
    K = h_mat.shape[0]
    n_snvs = h_mat.shape[1]
    top_n = max(1, int(0.05 * n_snvs))
    print(f"NMF k={K}, n SNVs = {n_snvs}, top-5% threshold = {top_n}")

    # Pre-compute per-component top-5% SNV positions
    top_pos = {}
    for c in range(K):
        idx = np.argsort(h_mat[c])[-top_n:]
        top_pos[c] = snv_pos[idx]

    rows = []
    for label, pos in TAG_SNPS.items():
        in_set = pos in set(int(p) for p in snv_pos)
        if in_set:
            j = np.where(snv_pos == pos)[0][0]
            for c in range(K):
                w = h_mat[c, j]
                # rank: position of this SNV in descending sort of component c
                rank = int((h_mat[c] > w).sum())   # 0 = top
                pctile = 100.0 * (rank + 1) / n_snvs
                # Nearest top-5% SNV
                dist = np.min(np.abs(top_pos[c] - pos))
                rows.append({
                    "tag_snp": label, "tag_pos": pos,
                    "in_NMF_set": True,
                    "component": c, "H_weight": float(w),
                    "rank_in_component": rank,
                    "rank_percentile": pctile,
                    "in_top5pct": (rank < top_n),
                    "nearest_top5pct_dist_bp": int(dist),
                })
        else:
            # Not in NMF set — find nearest SNVs (within 2 kb) per component
            distances_to_set = np.abs(snv_pos - pos)
            j_near = np.argmin(distances_to_set)
            nearest_pos = int(snv_pos[j_near])
            nearest_dist = int(distances_to_set[j_near])
            for c in range(K):
                w_near = h_mat[c, j_near]
                rank = int((h_mat[c] > w_near).sum())
                pctile = 100.0 * (rank + 1) / n_snvs
                rows.append({
                    "tag_snp": label, "tag_pos": pos,
                    "in_NMF_set": False,
                    "component": c, "H_weight": float(w_near),
                    "rank_in_component": rank,
                    "rank_percentile": pctile,
                    "in_top5pct": (rank < top_n),
                    "nearest_NMF_pos": nearest_pos,
                    "nearest_NMF_dist_bp": nearest_dist,
                })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "tag_snp_component_map.csv", index=False)
    print(f"\nwrote {OUT/'tag_snp_component_map.csv'}")

    # Print pivoted
    print("\n=== H weight per component (raw NMF loading) ===")
    piv = df.pivot(index="tag_snp", columns="component", values="H_weight")
    piv.columns = [f"c{c}" for c in piv.columns]
    print(piv.round(3).to_string())

    print("\n=== Rank percentile per component (lower = stronger loading) ===")
    piv = df.pivot(index="tag_snp", columns="component", values="rank_percentile")
    piv.columns = [f"c{c}" for c in piv.columns]
    print(piv.round(1).to_string())

    print("\n=== Top-5% membership (1 = in top 5%) ===")
    piv = df.pivot(index="tag_snp", columns="component",
                    values="in_top5pct").astype(int)
    piv.columns = [f"c{c}" for c in piv.columns]
    print(piv.to_string())

    print("\n=== Nearest top-5% SNV distance (bp) ===")
    if "nearest_top5pct_dist_bp" in df.columns:
        piv = df.pivot(index="tag_snp", columns="component",
                        values="nearest_top5pct_dist_bp")
        piv.columns = [f"c{c}" for c in piv.columns]
        print(piv.to_string())

    # ===== Figure: H weight heatmap, with rank percentile annotation =====
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
    })
    weight_pivot = df.pivot(index="tag_snp", columns="component",
                              values="H_weight")
    rank_pivot = df.pivot(index="tag_snp", columns="component",
                            values="rank_percentile")

    # Order tags by genomic position
    order = sorted(weight_pivot.index, key=lambda x: TAG_SNPS[x])
    weight_pivot = weight_pivot.reindex(order)
    rank_pivot = rank_pivot.reindex(order)

    fig, axes = plt.subplots(1, 2, figsize=(13, 3.6))

    # Panel A: H weight heatmap (normalized per component to [0, 1])
    ax = axes[0]
    H_norm = weight_pivot.values / np.nanmax(h_mat, axis=1).reshape(1, -1)
    im = ax.imshow(H_norm, aspect="auto", cmap="Reds",
                    vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(K)); ax.set_xticklabels([f"c{c}" for c in range(K)])
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=8)
    for i in range(len(order)):
        for j in range(K):
            v = H_norm[i, j]
            r = rank_pivot.values[i, j]
            txt = f"{weight_pivot.values[i, j]:.2f}\n[{r:.0f}%]"
            ax.text(j, i, txt, ha="center", va="center",
                     fontsize=6.8, color="#fff" if v > 0.6 else "#222")
    fig.colorbar(im, ax=ax, fraction=0.04, label="H / max H per component")
    ax.set_title("A. NMF H loading of tag SNPs (normalized; rank % below)",
                  loc="left", fontsize=10)

    # Panel B: rank percentile (lower = stronger; top 5% = below 5)
    ax = axes[1]
    rmat = rank_pivot.values
    im = ax.imshow(rmat, aspect="auto", cmap="Blues_r",
                    vmin=0, vmax=50, interpolation="nearest")
    ax.set_xticks(range(K)); ax.set_xticklabels([f"c{c}" for c in range(K)])
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=8)
    for i in range(len(order)):
        for j in range(K):
            v = rmat[i, j]
            ax.text(j, i, f"{v:.1f}%", ha="center", va="center",
                     fontsize=7, color="#fff" if v < 15 else "#222")
    fig.colorbar(im, ax=ax, fraction=0.04, label="rank percentile (lower = stronger)")
    # Mark top-5% boundary
    ax.set_title("B. rank percentile per component (≤5% = top signature)",
                  loc="left", fontsize=10)

    fig.suptitle(
        "Where do known tag SNPs live in NMF component space?",
        fontsize=11, y=1.05,
    )
    fig.tight_layout()
    out = OUT / "figure_tag_snp_component"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
