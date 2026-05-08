"""AFR-specific Leiden cluster deep-dive.

The 26-pop UMAP+Leiden run identified 12 super-pop-private clusters
(>=80% from one super-pop), of which 5 are AFR-private:
  L0  (n=179, AFR=93%)
  L12 (n=69,  AFR=86%)
  L19 (n=55,  AFR=91%)
  L21 (n=47,  AFR=87%)
  L31 (n=25,  AFR=92%)

Questions:
  1. Are these 5 AFR clusters the same "AFR background" or distinct lineages?
  2. Which NMF components load most heavily on each AFR cluster?
  3. Which AFR sub-pops (YRI/LWK/GWD/MSL/ESN/ASW/ACB) contribute to each?
  4. Spatial signature: which chromosomal regions distinguish them?

Output:
  results/afr_cluster_sub_pop_composition.csv
  results/afr_cluster_nmf_loadings.csv
  results/afr_cluster_signature_snvs.csv
  results/figure_afr_branches.{pdf,png}
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
VCF = REPO / "data/region_500kb.vcf.gz"
PANEL = Path("/home/yyamada1225/LVexome/demo_delta_mc/thousand_g/"
             "integrated_call_samples_v3.20130502.ALL.panel")
UMAP_PARQUET = REPO / "results/umap_leiden_26.parquet"
W_PARQUET = REPO / "results/nmf_W_26_k8.parquet"
H_PARQUET = REPO / "results/nmf_H_26_k8.parquet"
OUT = REPO / "results"

ANCHOR_POS = 31_366_595
AFR_POPS = ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"]
ALL_26 = {
    "AFR": AFR_POPS,
    "EUR": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "EAS": ["CHB", "JPT", "CHS", "CDX", "KHV"],
    "SAS": ["GIH", "PJL", "BEB", "STU", "ITU"],
    "AMR": ["MXL", "PUR", "CLM", "PEL"],
}
KEEP_POPS = [p for ps in ALL_26.values() for p in ps]
SUPER_OF = {p: s for s, ps in ALL_26.items() for p in ps}
MAF_FLOOR = 0.05


def load_matrix(path: Path):
    samples = None; positions, haps = [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]; continue
            f = line.rstrip("\n").split("\t")
            ref, alt = f[3], f[4]
            if (len(ref) != 1 or len(alt) != 1
                or alt.startswith("<") or "," in alt): continue
            pos = int(f[1])
            gts = f[9:]
            n = len(gts); hap = np.empty(2*n, dtype=np.int8)
            ok = True
            for i, g in enumerate(gts):
                if "|" not in g: ok = False; break
                a, b = g.split("|", 1)
                if a == "." or b == ".":
                    hap[2*i] = -1; hap[2*i+1] = -1
                else:
                    hap[2*i] = int(a); hap[2*i+1] = int(b)
            if not ok: continue
            positions.append(pos); haps.append(hap)
    return samples, np.array(positions), np.array(haps, dtype=np.int8)


def main():
    sys.stdout.reconfigure(line_buffering=True)
    umap_df = pd.read_parquet(UMAP_PARQUET)
    W = pd.read_parquet(W_PARQUET).values
    H = pd.read_parquet(H_PARQUET)
    snv_pos = np.array([int(c) for c in H.columns])
    h_mat = H.values

    # ---- AFR-private clusters from earlier output ----
    sup_x_leid = umap_df.groupby(["super_pop","leiden"]).size().unstack(fill_value=0)
    cluster_super_purity = sup_x_leid.div(sup_x_leid.sum(axis=0), axis=1)
    afr_private = [int(c) for c in sup_x_leid.columns
                   if cluster_super_purity[c].idxmax() == "AFR"
                   and cluster_super_purity[c].max() >= 0.80]
    afr_private.sort()
    sizes = sup_x_leid.sum(axis=0)
    print(f"AFR-private Leiden clusters: {afr_private}")
    for cl in afr_private:
        print(f"  L{cl}: n={int(sizes[cl])}, AFR purity={cluster_super_purity[cl]['AFR']:.0%}")

    # 1. Sub-pop composition
    print("\n=== AFR sub-pop composition (count) ===")
    afr_haps = umap_df[umap_df["super_pop"] == "AFR"]
    sub_table = (afr_haps[afr_haps["leiden"].isin(afr_private)]
                  .groupby(["leiden","pop"]).size().unstack(fill_value=0)
                  .reindex(columns=AFR_POPS, fill_value=0))
    print(sub_table.to_string())
    sub_freq = sub_table.div(sub_table.sum(axis=1), axis=0).round(3)
    print("\n=== AFR sub-pop composition (freq) ===")
    print(sub_freq.to_string())
    sub_table.to_csv(OUT / "afr_cluster_sub_pop_composition.csv")

    # 2. NMF component loadings per cluster
    print("\n=== Mean NMF component loading per AFR cluster ===")
    K = W.shape[1]
    rows = []
    for cl in afr_private:
        m = (umap_df["leiden"] == cl).values
        mean_w = W[m].mean(axis=0)
        rows.append({"cluster": f"L{cl}", "n_hap": int(m.sum()),
                      **{f"c{c}": float(mean_w[c]) for c in range(K)}})
    nmf_df = pd.DataFrame(rows)
    print(nmf_df.to_string(index=False, float_format="%.3f"))
    nmf_df.to_csv(OUT / "afr_cluster_nmf_loadings.csv", index=False)

    # 3. Cluster-distinguishing SNVs
    # Load full carrier matrix (re-derive to align with umap_df rows)
    print("\nloading full hap matrix...")
    samples, positions, haps = load_matrix(VCF)
    panel = pd.read_csv(PANEL, sep="\t")
    pop_of = dict(zip(panel["sample"], panel["pop"]))
    super_of = dict(zip(panel["sample"], panel["super_pop"]))
    hap_pop = np.empty(haps.shape[1], dtype=object)
    for i, s in enumerate(samples):
        hap_pop[2*i] = pop_of.get(s, "UNK"); hap_pop[2*i+1] = pop_of.get(s, "UNK")
    arow = np.where(positions == ANCHOR_POS)[0][0]
    is_carrier = (haps[arow] == 1) & np.isin(hap_pop, KEEP_POPS)
    cidx = np.where(is_carrier)[0]
    sub_haps = haps[:, cidx]
    mask = np.ones(len(positions), bool); mask[arow] = False
    sub_haps = sub_haps[mask]; sub_pos = positions[mask]
    valid = (sub_haps >= 0).all(axis=1)
    sub_haps = sub_haps[valid]; sub_pos = sub_pos[valid]
    af = sub_haps.mean(axis=1)
    info = (af >= MAF_FLOOR) & (af <= 1 - MAF_FLOOR)
    M_pos = sub_pos[info]
    M = sub_haps[info].T   # shape: 2111 x 7116, hap x SNV

    if M.shape[0] != len(umap_df):
        raise SystemExit("hap count mismatch with umap_df")

    # For each AFR cluster, find SNVs with maximally different freq vs other carriers
    print("\nfinding cluster-discriminating SNVs per AFR cluster (top 30)...")
    sig_rows = []
    for cl in afr_private:
        m_in = (umap_df["leiden"] == cl).values
        m_out = ~m_in
        af_in = M[m_in].mean(axis=0)
        af_out = M[m_out].mean(axis=0)
        diff = af_in - af_out
        top_idx = np.argsort(np.abs(diff))[-30:][::-1]
        for j in top_idx:
            sig_rows.append({
                "cluster": f"L{cl}", "partner_pos": int(M_pos[j]),
                "af_in": float(af_in[j]), "af_out": float(af_out[j]),
                "delta": float(diff[j]),
            })
    sig_df = pd.DataFrame(sig_rows)
    sig_df.to_csv(OUT / "afr_cluster_signature_snvs.csv", index=False)
    print(f"  wrote {OUT/'afr_cluster_signature_snvs.csv'}")

    # ---- Figures ----
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })

    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.35, wspace=0.3)

    # A: AFR sub-pop composition heatmap
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(sub_freq.values, aspect="auto", cmap="Blues",
                    vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(len(AFR_POPS))); ax.set_xticklabels(AFR_POPS,
                                                              rotation=30, ha="right")
    ax.set_yticks(range(len(sub_freq))); ax.set_yticklabels(sub_freq.index)
    for i in range(len(sub_freq)):
        for j in range(len(AFR_POPS)):
            v = sub_freq.values[i, j]
            if v > 0:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                         fontsize=7, color="#fff" if v > 0.5 else "#111")
    fig.colorbar(im, ax=ax, fraction=0.04)
    ax.set_title("A. AFR sub-pop composition", loc="left")

    # B: NMF component loading heatmap
    ax = fig.add_subplot(gs[0, 1])
    nmf_mat = nmf_df.iloc[:, 2:].values
    im = ax.imshow(nmf_mat, aspect="auto", cmap="viridis", interpolation="nearest")
    ax.set_xticks(range(K)); ax.set_xticklabels([f"c{c}" for c in range(K)])
    ax.set_yticks(range(len(nmf_df))); ax.set_yticklabels(nmf_df["cluster"])
    fig.colorbar(im, ax=ax, fraction=0.04)
    ax.set_title("B. Mean NMF component loading", loc="left")

    # C: Spatial signature for each AFR cluster
    ax = fig.add_subplot(gs[0, 2])
    afr_palette = plt.cm.Reds(np.linspace(0.4, 0.95, len(afr_private)))
    for i, cl in enumerate(afr_private):
        m_in = (umap_df["leiden"] == cl).values
        m_out = ~m_in
        af_in = M[m_in].mean(axis=0)
        af_out = M[m_out].mean(axis=0)
        diff = af_in - af_out
        ax.scatter(M_pos / 1e6, diff, s=2, color=afr_palette[i],
                    alpha=0.5, edgecolor="none", label=f"L{cl}")
    ax.axhline(0, color="#222", lw=0.5)
    for lpos, ln in [(31_366_595,"rs2596542"),(31_435_833,"rs2244546"),
                       (31_475_084,"rs11509487")]:
        ax.axvline(lpos/1e6, color="#222", lw=0.4, ls="--", alpha=0.4)
    ax.set_xlabel("chromosomal position (Mb, GRCh37)")
    ax.set_ylabel("AF (in cluster) − AF (other carriers)")
    ax.set_title("C. Cluster-discriminating SNVs", loc="left")
    ax.legend(fontsize=7, loc="upper right")

    # D-G: per-cluster spatial profile (top 4 clusters)
    for i, cl in enumerate(afr_private[:4]):
        ax = fig.add_subplot(gs[1, i] if i < 3 else gs[1, 2])
        m_in = (umap_df["leiden"] == cl).values
        m_out = ~m_in
        af_in = M[m_in].mean(axis=0)
        af_out = M[m_out].mean(axis=0)
        diff = af_in - af_out
        ax.scatter(M_pos/1e6, diff, s=3, color=afr_palette[i],
                    alpha=0.55, edgecolor="none")
        ax.axhline(0, color="#222", lw=0.4)
        for lpos, ln in [(31_366_595,"anchor"),(31_435_833,"rs2244546"),
                          (31_475_084,"rs11509487 (MICB)")]:
            ax.axvline(lpos/1e6, color="#222", lw=0.4, ls="--", alpha=0.4)
        ax.set_title(f"L{cl} (n={int(m_in.sum())})", loc="left", fontsize=9)
        ax.set_ylim(-0.6, 0.8)
        ax.set_xlabel("chr6 (Mb)")
        if i == 0: ax.set_ylabel("AF in cluster − AF out")
        if i == 3: break

    fig.suptitle(
        "AFR-specific anchor-carrier branches (rs2596542-T, 26-pop k=8 NMF)",
        fontsize=11, y=1.005,
    )
    out = OUT / "figure_afr_branches"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
