"""NMF k-sensitivity: do MICA-stable and HLA-B-variable axes emerge
consistently across k = {4, 5, 6, 7, 8, 9, 10, 12}?

For each k, we:
  1. fit NMF on the same hap × SNV matrix (rs2596542-T carriers, 26 pops);
  2. for each component, take top-5% SNVs;
  3. compute fold enrichment for: r-flip, C-flip, class-flip,
     and target-gene eQTL (MICA, MICB, HLA-B, HLA-C);
  4. auto-label component: MICA-stable / HLA-B-variable /
     HLA-B-stable / MICB / mixed.

Output:
  results/k_sensitivity_components.csv
  results/figure_k_sensitivity.{pdf,png}
"""

from __future__ import annotations

import gzip
import sys
import time
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from sklearn.decomposition import NMF

REPO = Path(__file__).resolve().parent
VCF = REPO / "data/region_500kb.vcf.gz"
PANEL = Path("/home/yyamada1225/LVexome/demo_delta_mc/thousand_g/"
             "integrated_call_samples_v3.20130502.ALL.panel")
EQTL_CSV = REPO / "results/eqtl_targets_per_gene.csv"
REGION_CSV = REPO / "results/region_per_partner_26.csv"
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
SUPER_POPS = list(ALL_26.keys())
MAF_FLOOR = 0.05
NMF_TOP_FRAC = 0.05
KS = [4, 5, 6, 7, 8, 9, 10, 12]


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
    return (2 * np.sum(np.arange(1, n+1) * x) / (n * x.sum())) - (n+1)/n


def label_component(row):
    """Auto-label based on enrichment profile."""
    mica  = row["MICA_fold"]
    hlab  = row["HLA-B_fold"]
    hlac  = row["HLA-C_fold"]
    micb  = row["MICB_fold"]
    rfl   = row["r_flip_fold"]
    cfl   = row["C_flip_fold"]
    clf   = row["class_flip_fold"]

    # MICA-stable: MICA >= 2, r/C flip <= 0.9, class flip > 1
    if mica >= 2.0 and (rfl <= 0.9) and (cfl <= 0.9):
        return "MICA-stable"
    # HLA-B-variable: HLA-B/C >= 1.5 + r/C flip >= 1.15, MICA <= 0.5
    if (hlab >= 1.5 or hlac >= 1.5) and (rfl >= 1.15 or cfl >= 1.15) and mica <= 0.5:
        return "HLA-B-variable"
    # HLA-B-stable: HLA-B >= 2 but r/C flip not enriched
    if hlab >= 2.0 and rfl < 1.15 and cfl < 1.15:
        return "HLA-B-stable"
    # MICB-broadcast: MICB >= 1.5, mild flip
    if micb >= 1.5 and rfl >= 1.10:
        return "MICB-broadcast"
    if mica < 0.3 and hlab < 1.0 and hlac < 1.0:
        return "low-signal"
    return "mixed"


def main():
    sys.stdout.reconfigure(line_buffering=True)
    t0 = time.time()

    # Load matrix once
    print("loading VCF...")
    samples, positions, haps = load_matrix(VCF)
    panel = pd.read_csv(PANEL, sep="\t")
    pop_of = dict(zip(panel["sample"], panel["pop"]))
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
    M = sub_haps[info].T.astype(np.float32)
    print(f"  matrix {M.shape}, {(time.time()-t0):.0f}s")

    # Load eQTL/flip sets
    eqtl = pd.read_csv(EQTL_CSV)
    eqtl_sets = {g: set(eqtl[eqtl["gene"] == g]["b37_pos"].astype(int))
                 for g in eqtl["gene"].unique()}
    region = pd.read_csv(REGION_CSV)
    universe = set(int(p) for p in M_pos)
    cflip_set = set(region.loc[region["any_C_flip"] == True,
                                  "partner_pos"].astype(int)) & universe
    rmat = region[[f"median_r_{s}" for s in SUPER_POPS]].astype(float).values
    nrf = np.zeros(len(region), dtype=int)
    for i, j in combinations(range(len(SUPER_POPS)), 2):
        ri = rmat[:, i]; rj = rmat[:, j]
        v = ~(np.isnan(ri) | np.isnan(rj))
        nrf += (v & ((np.sign(ri) * np.sign(rj)) < 0)).astype(int)
    rflip_set = set(region.loc[nrf >= 1, "partner_pos"].astype(int)) & universe
    classflip_set = set(region.loc[region["any_class_flip"] == True,
                                      "partner_pos"].astype(int)) & universe
    for g in eqtl_sets: eqtl_sets[g] &= universe

    print(f"  test sets: HLA-B={len(eqtl_sets.get('HLA-B', set()))}, "
           f"MICA={len(eqtl_sets.get('MICA', set()))}, "
           f"r-flip={len(rflip_set)}, C-flip={len(cflip_set)}, "
           f"class-flip={len(classflip_set)}")

    # ===== NMF for each k =====
    rows = []
    for k in KS:
        print(f"\nNMF k={k} ...")
        nmf = NMF(n_components=k, init="nndsvdar", random_state=0,
                   max_iter=400, tol=1e-4)
        W = nmf.fit_transform(M); H = nmf.components_
        n_snvs = H.shape[1]
        top_n = max(1, int(NMF_TOP_FRAC * n_snvs))
        for c in range(k):
            h = H[c]
            top_idx = np.argsort(h)[-top_n:]
            top_set = set(int(p) for p in M_pos[top_idx])
            row = {"k": k, "component": c, "n_top": top_n,
                    "hap_gini": gini(W[:, c]),
                    "snv_gini": gini(h)}
            for label, ts in [("r_flip", rflip_set), ("C_flip", cflip_set),
                                ("class_flip", classflip_set)]:
                K_succ = len(ts); N = len(universe); n_d = top_n
                x = len(top_set & ts)
                if K_succ == 0:
                    fold = float("nan"); p = float("nan")
                else:
                    fold = x / (n_d * K_succ / N) if K_succ > 0 else float("nan")
                    p = float(hypergeom.sf(x - 1, N, K_succ, n_d))
                row[f"{label}_fold"] = fold; row[f"{label}_p"] = p
            for gene in ["MICA", "MICB", "HLA-B", "HLA-C"]:
                ts = eqtl_sets.get(gene, set())
                K_succ = len(ts); N = len(universe); n_d = top_n
                x = len(top_set & ts)
                if K_succ == 0:
                    fold = float("nan"); p = float("nan")
                else:
                    fold = x / (n_d * K_succ / N)
                    p = float(hypergeom.sf(x - 1, N, K_succ, n_d))
                row[f"{gene}_fold"] = fold; row[f"{gene}_p"] = p
            row["label"] = label_component(row)
            rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "k_sensitivity_components.csv", index=False)
    print(f"\nwrote {OUT/'k_sensitivity_components.csv'}")

    # Print compact summary
    print("\n=== k × component label ===")
    pcounts = df.groupby(["k", "label"]).size().unstack(fill_value=0)
    print(pcounts.to_string())

    print("\n=== Components per k (HLA-B-variable + MICA-stable highlighted) ===")
    cols = ["k","component","label","MICA_fold","HLA-B_fold","HLA-C_fold",
             "MICB_fold","r_flip_fold","C_flip_fold","class_flip_fold"]
    print(df[cols].to_string(index=False, formatters={
        "MICA_fold":"{:.2f}".format, "HLA-B_fold":"{:.2f}".format,
        "HLA-C_fold":"{:.2f}".format, "MICB_fold":"{:.2f}".format,
        "r_flip_fold":"{:.2f}".format, "C_flip_fold":"{:.2f}".format,
        "class_flip_fold":"{:.2f}".format,
    }))

    # ===== Figure: per-k, scatter MICA fold vs HLA-B fold colored by C-flip fold =====
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharex=True, sharey=True)
    for ax, k in zip(axes.flat, KS):
        sub = df[df["k"] == k]
        sc = ax.scatter(sub["MICA_fold"], sub["HLA-B_fold"],
                         c=sub["C_flip_fold"], cmap="RdBu_r", vmin=0.6, vmax=1.4,
                         s=120, edgecolor="#222", linewidth=0.5)
        for _, r in sub.iterrows():
            label_short = {"MICA-stable":"M", "HLA-B-variable":"v",
                            "HLA-B-stable":"s", "MICB-broadcast":"b",
                            "mixed":"x", "low-signal":"."}.get(r["label"], "?")
            ax.annotate(f"c{r['component']}{label_short}",
                          (r["MICA_fold"], r["HLA-B_fold"]),
                          fontsize=7, ha="center", va="center")
        ax.axhline(1.0, color="#888", lw=0.4, ls="--")
        ax.axvline(1.0, color="#888", lw=0.4, ls="--")
        ax.set_title(f"k = {k}", loc="left")
        ax.set_xlabel("MICA eQTL fold")
        ax.set_ylabel("HLA-B eQTL fold")
        ax.grid(linestyle=":", alpha=0.3)
    cbar = fig.colorbar(sc, ax=axes, fraction=0.012, pad=0.01,
                         label="C-flip fold (red = enriched)")
    fig.suptitle(
        "k sensitivity — does the MICA-stable / HLA-B-variable axis dichotomy persist?\n"
        "(label codes: M=MICA-stable, v=HLA-B-variable, s=HLA-B-stable, b=MICB-broadcast, x=mixed)",
        fontsize=11, y=1.005,
    )
    out = OUT / "figure_k_sensitivity"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
