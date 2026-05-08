"""26-pop UMAP + Leiden + NMF + multi-super-pop C-flip enrichment.

Same pipeline as decompose_branches.py but extended to all 26 1000G
sub-populations (AFR + EUR + EAS + SAS + AMR).

C-flip definition (multi-super-pop):
  any_C_flip = at least one of 10 super-pop pairs has opposite-sign
  median(C). See analyze_region_MC_26.py.

Outputs:
  results/umap_leiden_26.parquet
  results/nmf_W_26_k{K}.parquet
  results/nmf_H_26_k{K}.parquet
  results/nmf_component_summary_26.csv
  results/cflip_enrichment_26.csv
  results/figure_umap_leiden_26.{pdf,png}
  results/figure_nmf_components_26.{pdf,png}
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
REGION_PARTNER_CSV = REPO / "results/region_per_partner_26.csv"
OUT = REPO / "results"

ANCHOR_POS = 31_366_595
ALL_26 = {
    "AFR": ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"],
    "EUR": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "EAS": ["CHB", "JPT", "CHS", "CDX", "KHV"],
    "SAS": ["GIH", "PJL", "BEB", "STU", "ITU"],
    "AMR": ["MXL", "PUR", "CLM", "PEL"],
}
KEEP_POPS = [p for ps in ALL_26.values() for p in ps]
SUPER_OF = {p: s for s, ps in ALL_26.items() for p in ps}
SUPER_POPS = list(ALL_26.keys())
MAF_FLOOR = 0.05
NMF_KS = [6, 8, 10]
LEIDEN_RES = 1.0
NMF_TOP_FRAC = 0.05


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


def gini(x):
    x = np.asarray(x, dtype=float)
    if x.sum() == 0: return 0.0
    x = np.sort(x); n = len(x)
    return (2 * np.sum(np.arange(1, n + 1) * x) / (n * x.sum())) - (n + 1) / n


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
    M = sub_haps[info].T.astype(np.float32)
    print(f"  matrix: {M.shape}  ({(time.time()-t0):.0f}s)")
    print(f"  carriers per super-pop:")
    for sup in SUPER_POPS:
        n = (sub_super == sup).sum()
        print(f"    {sup}: {n}")

    # ===== UMAP + Leiden =====
    print("\nUMAP...")
    reducer = umap.UMAP(n_neighbors=20, min_dist=0.05, metric="hamming",
                         random_state=0, n_components=2)
    umap_xy = reducer.fit_transform(M)
    print(f"  done ({(time.time()-t0):.0f}s)")

    print("Leiden...")
    nn = NearestNeighbors(n_neighbors=20, metric="hamming"); nn.fit(M)
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
    sizes = np.bincount(leiden_label).tolist()
    print(f"  {n_leiden} clusters, sizes {sizes[:10]}{'...' if len(sizes)>10 else ''} "
           f"({(time.time()-t0):.0f}s)")

    umap_df = pd.DataFrame({
        "pop": sub_pop, "super_pop": sub_super,
        "UMAP1": umap_xy[:, 0], "UMAP2": umap_xy[:, 1],
        "leiden": leiden_label,
    })
    umap_df.to_parquet(OUT / "umap_leiden_26.parquet")

    # Super-pop x leiden table
    print("\nLeiden cluster x super-pop (top clusters):")
    sup_x_leid = umap_df.groupby(["super_pop","leiden"]).size().unstack(fill_value=0)
    sup_x_leid_freq = sup_x_leid.div(sup_x_leid.sum(axis=1), axis=0)
    print(sup_x_leid_freq.round(3).to_string())

    # Identify pop-private clusters: those with >= 80% from one super-pop
    cluster_super_purity = sup_x_leid.div(sup_x_leid.sum(axis=0), axis=1)
    private = cluster_super_purity.max(axis=0) >= 0.8
    print(f"\n  super-pop-private clusters (>=80% from single super-pop): "
           f"{private.sum()} / {n_leiden}")
    private_ids = sup_x_leid.columns[private].tolist()
    for cl in private_ids[:10]:
        purity = cluster_super_purity[cl]
        dom = purity.idxmax()
        n = sup_x_leid[cl].sum()
        print(f"    L{cl} (n={n}): {dom} = {purity[dom]:.0%}")

    # ===== NMF =====
    nmf_results = {}
    for k in NMF_KS:
        print(f"\nNMF k={k}...")
        nmf = NMF(n_components=k, init="nndsvdar", random_state=0,
                  max_iter=400, tol=1e-4)
        W = nmf.fit_transform(M); H = nmf.components_
        rows = []
        for c in range(k):
            w = W[:, c]; h = H[c, :]
            rows.append({
                "k": k, "component": c,
                "hap_loading_gini": gini(w),
                "snv_signature_gini": gini(h),
                "n_top5pct_snvs": int(0.05 * len(h)),
            })
        nmf_results[k] = {"W": W, "H": H, "comp": pd.DataFrame(rows),
                            "err": nmf.reconstruction_err_}
        print(f"  reconstruction_err = {nmf.reconstruction_err_:.3f}")

    K_USE = 8
    pd.DataFrame(nmf_results[K_USE]["W"],
                  columns=[f"c{i}" for i in range(K_USE)]
                  ).to_parquet(OUT / f"nmf_W_26_k{K_USE}.parquet")
    pd.DataFrame(nmf_results[K_USE]["H"],
                  columns=[str(int(p)) for p in M_pos]
                  ).to_parquet(OUT / f"nmf_H_26_k{K_USE}.parquet")
    pd.concat([nmf_results[k]["comp"] for k in NMF_KS], ignore_index=True
               ).to_csv(OUT / "nmf_component_summary_26.csv", index=False)

    # ===== Multi-super-pop C-flip enrichment =====
    print("\nC-flip enrichment test (multi-super-pop)...")
    region = pd.read_csv(REGION_PARTNER_CSV)
    cflip_pos = set(region.loc[region["any_C_flip"] == True, "partner_pos"].astype(int))
    print(f"  total any-pair C-flip partners: {len(cflip_pos)}")

    nmf_pos_set = set(M_pos.astype(int))
    universe = nmf_pos_set & set(region["partner_pos"].astype(int))
    cflip_in_univ = cflip_pos & universe
    print(f"  universe (NMF ∩ region): {len(universe)}")
    print(f"  C-flip in universe: {len(cflip_in_univ)} "
           f"({len(cflip_in_univ)/len(universe)*100:.1f}%)")

    # Also test "high-flip" partners (>=4 pairs flip = at least 1 super-pop dissents)
    highflip_pos = set(region.loc[region["n_super_C_flips"] >= 4,
                                    "partner_pos"].astype(int))
    highflip_in_univ = highflip_pos & universe
    print(f"  >=4-pair flip partners: {len(highflip_pos)}, "
           f"in universe {len(highflip_in_univ)}")

    rows = []
    for k in NMF_KS:
        H = nmf_results[k]["H"]; W = nmf_results[k]["W"]
        for c in range(k):
            h = H[c, :]
            top_n = max(1, int(NMF_TOP_FRAC * len(h)))
            top_idx = np.argsort(h)[-top_n:]
            top_pos = set(M_pos[top_idx].astype(int))
            top_in_univ = top_pos & universe
            x_any = len(top_in_univ & cflip_pos)
            x_high = len(top_in_univ & highflip_pos)
            N = len(universe)
            n_draws = len(top_in_univ)
            K_any = len(cflip_in_univ); K_hi = len(highflip_in_univ)
            p_any = float(hypergeom.sf(x_any - 1, N, K_any, n_draws)) if n_draws else float("nan")
            p_hi  = float(hypergeom.sf(x_high - 1, N, K_hi, n_draws)) if n_draws else float("nan")
            exp_any = n_draws * K_any / N if N > 0 else float("nan")
            exp_hi  = n_draws * K_hi / N if N > 0 else float("nan")
            rows.append({
                "k": k, "component": c,
                "n_top": n_draws,
                "x_any_flip": x_any, "exp_any": exp_any,
                "fold_any": x_any / exp_any if exp_any > 0 else float("nan"),
                "p_any": p_any,
                "x_highflip": x_high, "exp_high": exp_hi,
                "fold_high": x_high / exp_hi if exp_hi > 0 else float("nan"),
                "p_high": p_hi,
                "hap_gini": gini(W[:, c]),
                "snv_gini": gini(h),
            })
    enrich = pd.DataFrame(rows)
    enrich.to_csv(OUT / "cflip_enrichment_26.csv", index=False)
    print()
    print(enrich.to_string(index=False, formatters={
        "exp_any": "{:.1f}".format,
        "fold_any": "{:.2f}".format,
        "p_any": "{:.2e}".format,
        "exp_high": "{:.1f}".format,
        "fold_high": "{:.2f}".format,
        "p_high": "{:.2e}".format,
        "hap_gini": "{:.3f}".format,
        "snv_gini": "{:.3f}".format,
    }))

    # ===== Figures =====
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
    })

    super_palette = {"AFR":"#b04545","EUR":"#3a8c5a","EAS":"#5a4fcf",
                      "SAS":"#b5722d","AMR":"#7e3da7"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    ax = axes[0]
    for sup in SUPER_POPS:
        m = umap_df["super_pop"] == sup
        ax.scatter(umap_df.loc[m, "UMAP1"], umap_df.loc[m, "UMAP2"],
                    s=8, color=super_palette[sup], alpha=0.55,
                    edgecolor="none", label=f"{sup} (n={m.sum()})")
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_title("A. UMAP — colored by super-pop", loc="left")
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(linestyle=":", alpha=0.3)

    ax = axes[1]
    leiden_palette = plt.cm.tab20(np.linspace(0, 1, max(20, n_leiden)))
    for cl in sorted(set(leiden_label)):
        m = umap_df["leiden"] == cl
        if m.sum() < 8: continue
        ax.scatter(umap_df.loc[m, "UMAP1"], umap_df.loc[m, "UMAP2"],
                    s=8, color=leiden_palette[cl % 20], alpha=0.7,
                    edgecolor="none")
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_title(f"B. Leiden clusters (n={n_leiden})", loc="left")
    ax.grid(linestyle=":", alpha=0.3)
    fig.suptitle(
        f"rs2596542-T anchor-carrier haplotypes — UMAP + Leiden, all 26 1000G pops "
        f"(n = {M.shape[0]}, {M.shape[1]} SNVs)",
        fontsize=10.5, y=1.02,
    )
    fig.tight_layout()
    out = OUT / "figure_umap_leiden_26"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)

    # NMF spatial figure
    H = nmf_results[K_USE]["H"]; W = nmf_results[K_USE]["W"]
    fig, axes = plt.subplots(K_USE, 1, figsize=(13, 1.2 * K_USE), sharex=True)
    for c in range(K_USE):
        h = H[c, :]
        top_n = max(1, int(NMF_TOP_FRAC * len(h)))
        top_idx = set(np.argsort(h)[-top_n:])
        is_top = np.array([i in top_idx for i in range(len(h))])
        col = np.where(is_top, "#d65454", "#888888")
        axes[c].scatter(M_pos/1e6, h, s=4, c=col, alpha=0.6, edgecolor="none")
        gw = gini(W[:, c]); gh = gini(h)
        axes[c].set_ylabel(f"c{c}\nGw={gw:.2f}\nGh={gh:.2f}", fontsize=7.5)
        axes[c].grid(linestyle=":", alpha=0.3)
        for lpos, ln in [(31_366_595,"rs2596542"),(31_435_833,"rs2244546"),
                          (31_475_084,"rs11509487")]:
            axes[c].axvline(lpos/1e6, color="#222", lw=0.4, ls="--", alpha=0.5)
    axes[-1].set_xlabel("chromosomal position (Mb, GRCh37)")
    fig.suptitle(f"NMF k={K_USE} — SNV signatures (top 5% red), 26-pop carriers",
                  fontsize=10.5, y=1.005)
    fig.tight_layout()
    out = OUT / "figure_nmf_components_26"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"wrote {out}")
    plt.close(fig)

    print(f"\nTotal: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
