"""Cluster anchor-carrier haplotypes around rs2596542.

Hypothesis:
  rs2596542 minor allele (T, GRCh37 forward) tags a UNION of multiple
  haplotype branches. Branch composition differs across populations,
  driving the C-flip pattern observed in the region-wide scan.

Procedure:
  1. Load phased haplotypes for chr6:31,116,595-31,616,595 (±250 kb).
  2. Subset to haplotypes carrying rs2596542-T (the risk allele in the
     Kumar/Huang convention; forward-strand T).
  3. Build binary haplotype × SNV matrix; restrict SNVs to MAF >= 0.05
     within anchor carriers (informative variants only).
  4. Hierarchical clustering (Hamming + average linkage); cut at k=2,3,4,5.
  5. PCA for visual sanity check.
  6. Per-population branch composition.

Output:
  results/anchor_haplotype_branches.parquet (haplotype-level labels)
  results/branch_pop_composition.csv (k=2..5 frequencies)
  results/figure_branch_pca.{pdf,png}
  results/figure_branch_composition.{pdf,png}
"""

from __future__ import annotations

import gzip
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA

REPO = Path(__file__).resolve().parent
VCF = REPO / "data/region_500kb.vcf.gz"
PANEL = Path("/home/yyamada1225/LVexome/demo_delta_mc/thousand_g/"
             "integrated_call_samples_v3.20130502.ALL.panel")
OUT = REPO / "results"
OUT.mkdir(exist_ok=True)

ANCHOR_POS = 31_366_595  # rs2596542
EAS_POPS = ["CHB", "JPT", "CHS", "CDX", "KHV"]
EUR_POPS = ["CEU", "TSI", "FIN", "GBR", "IBS"]
KEEP_POPS = EAS_POPS + EUR_POPS
MAF_FLOOR = 0.05  # within anchor carriers


def load_haplotypes(path: Path):
    samples = None
    positions = []
    haplotypes = []  # list of np.int8 arrays (one per SNV)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            ref, alt = f[3], f[4]
            if (len(ref) != 1 or len(alt) != 1
                or alt.startswith("<") or "," in alt):
                continue
            pos = int(f[1])
            gts = f[9:]
            n = len(gts)
            hap = np.empty(2 * n, dtype=np.int8)
            ok = True
            for i, g in enumerate(gts):
                if "|" not in g:
                    ok = False
                    break
                a, b = g.split("|", 1)
                if a == "." or b == ".":
                    hap[2*i] = -1
                    hap[2*i+1] = -1
                else:
                    hap[2*i] = int(a)
                    hap[2*i+1] = int(b)
            if not ok:
                continue
            positions.append(pos)
            haplotypes.append(hap)
    return samples, np.array(positions), np.array(haplotypes, dtype=np.int8)


def main():
    sys.stdout.reconfigure(line_buffering=True)
    t0 = time.time()
    print("loading VCF...")
    samples, positions, haps = load_haplotypes(VCF)
    print(f"  {len(samples)} samples, {len(positions)} SNVs, "
           f"{haps.shape[1]} haplotypes  ({time.time()-t0:.0f}s)")

    panel = pd.read_csv(PANEL, sep="\t")
    pop_of = dict(zip(panel["sample"], panel["pop"]))
    super_of = dict(zip(panel["sample"], panel["super_pop"]))

    # Per-haplotype population
    hap_pop = np.empty(haps.shape[1], dtype=object)
    hap_super = np.empty(haps.shape[1], dtype=object)
    for i, s in enumerate(samples):
        hap_pop[2*i] = pop_of.get(s, "UNK")
        hap_pop[2*i+1] = pop_of.get(s, "UNK")
        hap_super[2*i] = super_of.get(s, "UNK")
        hap_super[2*i+1] = super_of.get(s, "UNK")

    # Find anchor row
    anchor_row = np.where(positions == ANCHOR_POS)[0]
    if len(anchor_row) == 0:
        raise SystemExit(f"anchor at {ANCHOR_POS} not in VCF")
    anchor_row = anchor_row[0]
    anchor_alleles = haps[anchor_row]

    # Anchor carriers: rs2596542-T (alt = 1 in dbSNP forward)
    is_carrier = (anchor_alleles == 1)
    # Restrict to EAS+EUR sub-pops (10 pops)
    in_pops = np.isin(hap_pop, KEEP_POPS)
    keep = is_carrier & in_pops
    n_anchor = keep.sum()
    print(f"  anchor carriers (rs2596542-T) in 10 pops: {n_anchor}")
    print(f"  per-pop counts:")
    for pop in KEEP_POPS:
        nk = ((hap_pop == pop) & keep).sum()
        print(f"    {pop}: {nk}")

    # Subset haplotypes
    sub_haps = haps[:, keep]
    sub_pop = hap_pop[keep]
    sub_super = hap_super[keep]

    # Drop the anchor row itself (no longer informative)
    mask_keep_row = np.ones(len(positions), dtype=bool)
    mask_keep_row[anchor_row] = False
    sub_haps = sub_haps[mask_keep_row]
    sub_positions = positions[mask_keep_row]

    # MAF filter within anchor carriers
    valid = (sub_haps >= 0).all(axis=1)
    sub_haps = sub_haps[valid]
    sub_positions = sub_positions[valid]
    af = sub_haps.mean(axis=1)
    informative = (af >= MAF_FLOOR) & (af <= 1 - MAF_FLOOR)
    M = sub_haps[informative].T  # shape: n_haplotypes x n_snv (informative)
    M_pos = sub_positions[informative]
    print(f"  informative SNVs (MAF>={MAF_FLOOR} in carriers): {M.shape[1]}")
    print(f"  matrix shape: {M.shape}")

    # Hierarchical clustering on Hamming distance
    print(f"  computing Hamming distance ({time.time()-t0:.0f}s)...")
    D = pdist(M, metric="hamming")
    print(f"    pdist done ({time.time()-t0:.0f}s)")
    Z = linkage(D, method="average")
    print(f"    linkage done ({time.time()-t0:.0f}s)")

    # Cut at k=2..5
    branch_labels = {}
    for k in [2, 3, 4, 5]:
        branch_labels[k] = fcluster(Z, t=k, criterion="maxclust")

    # PCA on the same matrix (centered)
    print(f"  PCA ({time.time()-t0:.0f}s)...")
    pca = PCA(n_components=2)
    XY = pca.fit_transform(M.astype(np.float32))
    print(f"    explained variance ratio: {pca.explained_variance_ratio_}")

    # Save haplotype-level table
    df = pd.DataFrame({
        "pop": sub_pop, "super_pop": sub_super,
        "PC1": XY[:, 0], "PC2": XY[:, 1],
        **{f"branch_k{k}": branch_labels[k] for k in [2, 3, 4, 5]}
    })
    df.to_parquet(OUT / "anchor_haplotype_branches.parquet")
    print(f"  wrote {OUT/'anchor_haplotype_branches.parquet'}  "
           f"({len(df)} haplotypes)")

    # Branch × pop composition (k=2..5)
    summary_rows = []
    for k in [2, 3, 4, 5]:
        for pop in KEEP_POPS:
            sub = df[df["pop"] == pop]
            n = len(sub)
            for br in range(1, k + 1):
                cnt = (sub[f"branch_k{k}"] == br).sum()
                summary_rows.append({
                    "k": k, "branch": br, "pop": pop,
                    "super_pop": "EAS" if pop in EAS_POPS else "EUR",
                    "n_carriers": n,
                    "n_in_branch": int(cnt),
                    "freq_in_branch": cnt / n if n > 0 else float("nan"),
                })
    comp = pd.DataFrame(summary_rows)
    comp.to_csv(OUT / "branch_pop_composition.csv", index=False)
    print(f"  wrote {OUT/'branch_pop_composition.csv'}")

    # ===== Figures =====
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
    })

    # Figure A: PCA scatter colored by pop and by branch (k=3)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    pop_palette = {
        "CHB": "#1f4e8c", "JPT": "#3070b8", "CHS": "#5089c8",
        "CDX": "#70a3d8", "KHV": "#90badd",
        "CEU": "#1f6c2e", "TSI": "#3a8c5a", "FIN": "#5aa874",
        "GBR": "#7ac28e", "IBS": "#9adba8",
    }
    ax = axes[0]
    for pop in KEEP_POPS:
        m = df["pop"] == pop
        ax.scatter(df.loc[m, "PC1"], df.loc[m, "PC2"],
                    s=12, color=pop_palette[pop], alpha=0.6,
                    edgecolor="none", label=pop)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title("A. anchor-carrier haplotypes (PCA), colored by pop", loc="left")
    ax.legend(fontsize=6.5, ncol=2, loc="upper right")
    ax.grid(linestyle=":", alpha=0.3)

    ax = axes[1]
    branch_palette = ["#d65454", "#3070b8", "#3a8c5a", "#b5722d", "#7e3da7"]
    k = 3
    for br in range(1, k + 1):
        m = df[f"branch_k{k}"] == br
        ax.scatter(df.loc[m, "PC1"], df.loc[m, "PC2"],
                    s=12, color=branch_palette[br - 1], alpha=0.6,
                    edgecolor="none", label=f"branch H{br}")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title(f"B. same haplotypes, colored by branch (k={k})", loc="left")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(linestyle=":", alpha=0.3)

    fig.suptitle(
        f"rs2596542-T anchor-carrier haplotypes "
        f"(n = {n_anchor} haplotypes, {M.shape[1]} informative SNVs, ±250 kb)",
        fontsize=10.5, y=1.02,
    )
    fig.tight_layout()
    out = OUT / "figure_branch_pca"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"  wrote {out}")
    plt.close(fig)

    # Figure B: branch composition per population, k=2..5
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.4), sharey=True)
    for ax, k in zip(axes, [2, 3, 4, 5]):
        wide = comp[comp["k"] == k].pivot(
            index="pop", columns="branch", values="freq_in_branch"
        ).reindex(KEEP_POPS)
        bottom = np.zeros(len(wide))
        x = np.arange(len(wide))
        for br in range(1, k + 1):
            vals = wide[br].values
            ax.bar(x, vals, bottom=bottom,
                    color=branch_palette[br - 1], edgecolor="#222",
                    linewidth=0.4, label=f"H{br}")
            bottom += vals
        ax.set_xticks(x)
        ax.set_xticklabels(wide.index, rotation=45, ha="right", fontsize=7.5)
        ax.set_title(f"k = {k}", fontsize=10, loc="left")
        if ax is axes[0]:
            ax.set_ylabel("branch frequency within rs2596542-T carriers")
        ax.set_ylim(0, 1)
        ax.axvline(4.5, color="#666", lw=0.8, ls="--")  # EAS|EUR separator
        ax.legend(fontsize=7, loc="upper right", ncol=1)
        ax.grid(axis="y", linestyle=":", alpha=0.3)
    fig.suptitle("Branch composition by population — rs2596542-T anchor",
                  fontsize=11, y=1.02)
    fig.tight_layout()
    out = OUT / "figure_branch_composition"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"  wrote {out}")
    plt.close(fig)

    print(f"\nTotal: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
