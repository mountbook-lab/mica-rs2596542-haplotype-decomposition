"""UMAP + Leiden + NMF decomposition of rs2596542-T anchor-carrier haplotypes.

Goal:
  Decompose anchor-carrier haplotype variation into:
    - branch-specific components (sparse hap loading, concentrated SNV signature)
    - backbone components (broad hap loading, spread SNV signature)
  and test whether C-flip partners (from region scan) preferentially fall
  into branch-specific signatures.

Method:
  1. Load 675 hap x 7,303 informative SNV matrix (built from VCF).
  2. UMAP 2D + Leiden clustering on kNN graph -- data-driven branch count.
  3. NMF (k = 4, 6, 8) on the binary matrix.
  4. For each NMF component:
        - hap-loading Gini (sparsity): high = branch-specific
        - SNV-signature Gini (concentration): high = focused signature
        - top 5 % SNVs by H-weight = component "tag set"
  5. Hypergeometric test: is the C-flip partner set
     (region_per_partner.csv with sign(C) flip EAS vs EUR)
     enriched in each component's tag set?

Outputs:
  results/umap_leiden.parquet
  results/nmf_W.parquet  (hap x component loadings)
  results/nmf_H.parquet  (component x SNV signatures)
  results/nmf_component_summary.csv
  results/cflip_enrichment.csv
  results/figure_umap_leiden.{pdf,png}
  results/figure_nmf_components.{pdf,png}
"""

from __future__ import annotations

import gzip
import sys
import time
from pathlib import Path

import igraph as ig
import leidenalg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
from scipy.stats import hypergeom
from sklearn.decomposition import NMF
from sklearn.neighbors import NearestNeighbors

REPO = Path(__file__).resolve().parent
VCF = REPO / "data/region_500kb.vcf.gz"
PANEL = Path("/home/yyamada1225/LVexome/demo_delta_mc/thousand_g/"
             "integrated_call_samples_v3.20130502.ALL.panel")
REGION_PARTNER_CSV = REPO / "results/region_per_partner.csv"
OUT = REPO / "results"

ANCHOR_POS = 31_366_595
EAS_POPS = ["CHB", "JPT", "CHS", "CDX", "KHV"]
EUR_POPS = ["CEU", "TSI", "FIN", "GBR", "IBS"]
KEEP_POPS = EAS_POPS + EUR_POPS
MAF_FLOOR = 0.05
NMF_KS = [4, 6, 8]
LEIDEN_RES = 1.0
NMF_TOP_FRAC = 0.05  # top 5% SNVs by H weight = component tag set


def load_matrix(path: Path):
    samples = None
    positions, haps = [], []
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


def gini(x):
    """Gini coefficient on a non-negative vector."""
    x = np.asarray(x, dtype=float)
    if x.sum() == 0: return 0.0
    x = np.sort(x); n = len(x)
    return (2 * np.sum((np.arange(1, n + 1)) * x) / (n * x.sum())) - (n + 1) / n


def main():
    sys.stdout.reconfigure(line_buffering=True)
    t0 = time.time()
    print("loading VCF...")
    samples, positions, haps = load_matrix(VCF)
    panel = pd.read_csv(PANEL, sep="\t")
    pop_of = dict(zip(panel["sample"], panel["pop"]))
    super_of = dict(zip(panel["sample"], panel["super_pop"]))
    hap_pop = np.empty(haps.shape[1], dtype=object)
    hap_super = np.empty(haps.shape[1], dtype=object)
    for i, s in enumerate(samples):
        hap_pop[2*i] = pop_of.get(s, "UNK"); hap_pop[2*i+1] = pop_of.get(s, "UNK")
        hap_super[2*i] = super_of.get(s, "UNK"); hap_super[2*i+1] = super_of.get(s, "UNK")

    arow = np.where(positions == ANCHOR_POS)[0][0]
    is_carrier = (haps[arow] == 1) & np.isin(hap_pop, KEEP_POPS)
    cidx = np.where(is_carrier)[0]
    sub_pop = hap_pop[cidx]; sub_super = hap_super[cidx]
    sub_haps = haps[:, cidx]
    mask = np.ones(len(positions), bool); mask[arow] = False
    sub_haps = sub_haps[mask]; sub_pos = positions[mask]
    valid = (sub_haps >= 0).all(axis=1)
    sub_haps = sub_haps[valid]; sub_pos = sub_pos[valid]
    af = sub_haps.mean(axis=1)
    info = (af >= MAF_FLOOR) & (af <= 1 - MAF_FLOOR)
    M_pos = sub_pos[info]
    M = sub_haps[info].T.astype(np.float32)   # (675, 7303) hap x SNV
    print(f"  matrix: {M.shape}, {(time.time()-t0):.0f}s")

    # ===== UMAP + Leiden =====
    print("\nUMAP...")
    reducer = umap.UMAP(n_neighbors=20, min_dist=0.05, metric="hamming",
                         random_state=0, n_components=2)
    umap_xy = reducer.fit_transform(M)
    print(f"  done ({(time.time()-t0):.0f}s)")

    print("Leiden...")
    nn = NearestNeighbors(n_neighbors=20, metric="hamming")
    nn.fit(M)
    knn = nn.kneighbors_graph(M, mode="connectivity")
    src, dst = knn.nonzero()
    edges = [(int(s), int(d)) for s, d in zip(src, dst) if s < d]
    g = ig.Graph(n=M.shape[0], edges=edges, directed=False)
    part = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=LEIDEN_RES, seed=0,
    )
    leiden_label = np.array(part.membership)
    n_leiden = len(set(leiden_label))
    print(f"  {n_leiden} clusters, sizes: "
           f"{np.bincount(leiden_label).tolist()} ({(time.time()-t0):.0f}s)")

    umap_df = pd.DataFrame({
        "pop": sub_pop, "super_pop": sub_super,
        "UMAP1": umap_xy[:, 0], "UMAP2": umap_xy[:, 1],
        "leiden": leiden_label,
    })
    umap_df.to_parquet(OUT / "umap_leiden.parquet")
    print(f"  wrote {OUT/'umap_leiden.parquet'}")

    # Per-leiden cluster pop composition
    print("\nLeiden cluster x population:")
    pop_x_leid = umap_df.groupby(["pop", "leiden"]).size().unstack(fill_value=0)
    pop_x_leid = pop_x_leid.reindex(KEEP_POPS)
    pop_x_leid_freq = pop_x_leid.div(pop_x_leid.sum(axis=1), axis=0)
    print(pop_x_leid_freq.round(2).to_string())

    # ===== NMF =====
    nmf_results = {}
    for k in NMF_KS:
        print(f"\nNMF k={k}...")
        nmf = NMF(n_components=k, init="nndsvdar", random_state=0,
                  max_iter=400, tol=1e-4)
        W = nmf.fit_transform(M); H = nmf.components_
        recon_err = nmf.reconstruction_err_
        # Per-component metrics
        rows = []
        for c in range(k):
            w = W[:, c]; h = H[c, :]
            rows.append({
                "k": k, "component": c,
                "hap_loading_gini": gini(w),
                "snv_signature_gini": gini(h),
                "hap_loading_top1pct_share":
                    np.sort(w)[-max(1, int(0.01*len(w))):].sum() / (w.sum() + 1e-9),
                "snv_top1pct_share":
                    np.sort(h)[-max(1, int(0.01*len(h))):].sum() / (h.sum() + 1e-9),
                "n_top5pct_snvs": int(0.05 * len(h)),
            })
        comp_df = pd.DataFrame(rows)
        nmf_results[k] = {"W": W, "H": H, "comp": comp_df, "err": recon_err}
        print(f"  reconstruction_err = {recon_err:.3f}")
        print(comp_df.to_string(index=False, formatters={
            "hap_loading_gini": "{:.3f}".format,
            "snv_signature_gini": "{:.3f}".format,
            "hap_loading_top1pct_share": "{:.3f}".format,
            "snv_top1pct_share": "{:.3f}".format,
        }))

    # Save W/H for chosen k (use k=6 by default)
    K_USE = 6
    pd.DataFrame(nmf_results[K_USE]["W"],
                  columns=[f"c{i}" for i in range(K_USE)]
                  ).to_parquet(OUT / f"nmf_W_k{K_USE}.parquet")
    pd.DataFrame(nmf_results[K_USE]["H"],
                  columns=[str(int(p)) for p in M_pos]
                  ).to_parquet(OUT / f"nmf_H_k{K_USE}.parquet")
    pd.concat([nmf_results[k]["comp"] for k in NMF_KS], ignore_index=True
               ).to_csv(OUT / "nmf_component_summary.csv", index=False)

    # ===== C-flip enrichment =====
    print("\nC-flip enrichment test...")
    region = pd.read_csv(REGION_PARTNER_CSV)
    region["c_flip"] = (
        ((region["median_C_EAS"] > 0) & (region["median_C_EUR"] < 0))
        | ((region["median_C_EAS"] < 0) & (region["median_C_EUR"] > 0))
    )
    cflip_pos = set(region.loc[region["c_flip"], "partner_pos"].astype(int))
    print(f"  total C-flip partners (region scan): {len(cflip_pos)}")

    # Use M_pos (NMF SNV positions) ∩ region partners as the universe
    nmf_pos_set = set(M_pos.astype(int))
    universe = nmf_pos_set & set(region["partner_pos"].astype(int))
    universe_arr = np.array(sorted(universe))
    cflip_in_univ = cflip_pos & universe
    print(f"  universe (in both NMF and region scan): {len(universe)}")
    print(f"  C-flip in universe: {len(cflip_in_univ)} "
          f"({len(cflip_in_univ)/len(universe)*100:.1f}%)")

    enrich_rows = []
    for k in NMF_KS:
        H = nmf_results[k]["H"]; W = nmf_results[k]["W"]
        for c in range(k):
            h = H[c, :]
            top_n = max(1, int(NMF_TOP_FRAC * len(h)))
            top_idx = np.argsort(h)[-top_n:]
            top_pos = set(M_pos[top_idx].astype(int))
            top_in_univ = top_pos & universe
            top_and_cflip = top_in_univ & cflip_pos
            # Hypergeometric test
            N = len(universe)
            K_succ = len(cflip_in_univ)
            n_draws = len(top_in_univ)
            x = len(top_and_cflip)
            if n_draws == 0 or N == 0:
                p_val = float("nan"); odds = float("nan")
            else:
                # P(X >= x)
                p_val = float(hypergeom.sf(x - 1, N, K_succ, n_draws))
                expected = n_draws * K_succ / N
                odds = x / expected if expected > 0 else float("nan")
            enrich_rows.append({
                "k": k, "component": c,
                "n_top_snvs_in_univ": n_draws,
                "n_cflip_in_top": x,
                "expected_cflip_in_top": n_draws * K_succ / N if N > 0 else float("nan"),
                "fold_enrichment": odds,
                "hyper_p": p_val,
                "hap_loading_gini": gini(W[:, c]),
                "snv_signature_gini": gini(h),
            })
    enrich = pd.DataFrame(enrich_rows)
    enrich.to_csv(OUT / "cflip_enrichment.csv", index=False)
    print()
    print(enrich.to_string(index=False, formatters={
        "expected_cflip_in_top": "{:.1f}".format,
        "fold_enrichment": "{:.2f}".format,
        "hyper_p": "{:.2e}".format,
        "hap_loading_gini": "{:.3f}".format,
        "snv_signature_gini": "{:.3f}".format,
    }))

    # ===== Figure: UMAP =====
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
    })
    pop_palette = {
        "CHB":"#1f4e8c","JPT":"#3070b8","CHS":"#5089c8","CDX":"#70a3d8","KHV":"#90badd",
        "CEU":"#1f6c2e","TSI":"#3a8c5a","FIN":"#5aa874","GBR":"#7ac28e","IBS":"#9adba8",
    }
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    ax = axes[0]
    for pop in KEEP_POPS:
        m = umap_df["pop"] == pop
        ax.scatter(umap_df.loc[m, "UMAP1"], umap_df.loc[m, "UMAP2"],
                    s=12, color=pop_palette[pop], alpha=0.7,
                    edgecolor="none", label=pop)
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_title("A. UMAP — colored by population", loc="left")
    ax.legend(fontsize=6.5, ncol=2, loc="upper right")
    ax.grid(linestyle=":", alpha=0.3)

    ax = axes[1]
    leiden_palette = plt.cm.tab10(np.linspace(0, 1, max(10, n_leiden)))
    for cl in sorted(set(leiden_label)):
        m = umap_df["leiden"] == cl
        ax.scatter(umap_df.loc[m, "UMAP1"], umap_df.loc[m, "UMAP2"],
                    s=12, color=leiden_palette[cl % 10], alpha=0.75,
                    edgecolor="none", label=f"L{cl} (n={m.sum()})")
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_title(f"B. Leiden clusters (k={n_leiden})", loc="left")
    ax.legend(fontsize=7, loc="upper right", ncol=1)
    ax.grid(linestyle=":", alpha=0.3)
    fig.suptitle(
        f"rs2596542-T anchor-carrier haplotypes — UMAP + Leiden "
        f"(n = {M.shape[0]}, {M.shape[1]} SNVs)",
        fontsize=10.5, y=1.02,
    )
    fig.tight_layout()
    out = OUT / "figure_umap_leiden"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)

    # ===== Figure: NMF components (k=K_USE) =====
    H = nmf_results[K_USE]["H"]; W = nmf_results[K_USE]["W"]
    fig, axes = plt.subplots(K_USE, 1, figsize=(13, 1.3 * K_USE), sharex=True)
    if K_USE == 1: axes = [axes]
    for c in range(K_USE):
        h = H[c, :]
        # mark top 5 % SNVs in red
        top_n = max(1, int(NMF_TOP_FRAC * len(h)))
        top_idx = set(np.argsort(h)[-top_n:])
        is_top = np.array([i in top_idx for i in range(len(h))])
        col = np.where(is_top, "#d65454", "#888888")
        axes[c].scatter(M_pos/1e6, h, s=4, c=col, alpha=0.6, edgecolor="none")
        gini_w = gini(W[:, c]); gini_h = gini(h)
        ftype = "branch" if (gini_w > 0.40 and gini_h > 0.55) else "backbone"
        axes[c].set_ylabel(f"c{c}\n({ftype}\nGw={gini_w:.2f}\nGh={gini_h:.2f})",
                            fontsize=7.5)
        axes[c].grid(linestyle=":", alpha=0.3)
        # landmarks
        for lpos, lname in [(31_366_595, "rs2596542"),
                              (31_435_833, "rs2244546"),
                              (31_475_084, "rs11509487")]:
            axes[c].axvline(lpos/1e6, color="#222", lw=0.4, ls="--", alpha=0.5)
    axes[-1].set_xlabel("chromosomal position (Mb, GRCh37)")
    fig.suptitle(
        f"NMF k={K_USE} components — SNV signatures (top 5% in red)",
        fontsize=10.5, y=1.005,
    )
    fig.tight_layout()
    out = OUT / "figure_nmf_components"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"wrote {out}")
    plt.close(fig)

    print(f"\nTotal: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
