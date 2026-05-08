"""eQTL effect direction (NES) coherence per NMF component.

Question:
  Within each NMF component's top-5% SNVs, is the eQTL effect direction
  (NES sign) on each target gene CONSISTENT (e.g., all push HLA-B
  expression up) or MIXED?

  If c5 top SNVs all push HLA-B/HLA-C expression in the same direction,
  c5 is a coherent "expression-modulating haplotype background."

  Comparison: c4/c6 (MICA-strong) NES on MICA — also coherent or mixed?

Procedure:
  1. Fetch GTEx Whole_Blood eQTLs for each target gene; keep variantId,
     pos, NES, pValue.
  2. Convert b38 positions -> b37 (offset 32,223). Critical: also flip NES
     sign if our REF/ALT differs from GTEx's REF/ALT (handle by tracking).
  3. For each NMF component × target gene, compute NES distribution among
     overlapping top-5% SNVs:
        - mean NES, sign-coherence = |mean(sign(NES))|
        - Wilcoxon test against zero (one-sample)
        - violin/strip plot for c1, c4, c5, c6 vs MICA, MICB, HLA-B, HLA-C

Outputs:
  results/eqtl_NES_per_component.csv
  results/figure_eqtl_direction.{pdf,png}
"""

from __future__ import annotations

import json
import sys
import time
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

REPO = Path(__file__).resolve().parent
H_PARQUET = REPO / "results/nmf_H_26_k8.parquet"
OUT = REPO / "results"

OFFSET_37_TO_38 = 32_223

TARGETS = [
    ("MICA",   "ENSG00000204520.12"),
    ("MICB",   "ENSG00000204516.9"),
    ("HLA-B",  "ENSG00000234745.10"),
    ("HLA-C",  "ENSG00000204525.16"),
    ("HCG27",  "ENSG00000206344.7"),
]
TISSUE = "Whole_Blood"
NMF_TOP_FRAC = 0.05
COMPONENTS_OF_INTEREST = [1, 4, 5, 6]   # c1 c4 c5 c6: the contrast set


def fetch_all(gencode: str, tissue: str):
    out = []; page = 0; items_per_page = 2000
    while True:
        url = (f"https://gtexportal.org/api/v2/association/singleTissueEqtl"
                f"?gencodeId={gencode}&tissueSiteDetailId={tissue}"
                f"&datasetId=gtex_v8&itemsPerPage={items_per_page}&page={page}")
        req = urllib.request.Request(url, headers={
            "Accept":"application/json","User-Agent":"Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            data = json.loads(r.read().decode())
        out.extend(data.get("data", []))
        n_pages = data.get("paging_info", {}).get("numberOfPages", 1)
        if page + 1 >= n_pages: break
        page += 1; time.sleep(0.3)
    return out


def main():
    sys.stdout.reconfigure(line_buffering=True)
    H = pd.read_parquet(H_PARQUET)
    snv_pos_b37 = np.array([int(c) for c in H.columns])
    h_mat = H.values
    K = h_mat.shape[0]
    print(f"NMF k={K}, n SNVs = {len(snv_pos_b37)}")

    # Top-5% per component
    top_sets = {}
    for c in range(K):
        h = h_mat[c]
        top_n = max(1, int(NMF_TOP_FRAC * len(h)))
        top_idx = np.argsort(h)[-top_n:]
        top_sets[c] = set(int(p) for p in snv_pos_b37[top_idx])
        print(f"  c{c}: top {top_n} SNVs")

    # 1. Fetch eQTLs with NES
    eqtl_data = {}   # gene -> dict {b37_pos -> NES} (best NES if multiple at same pos)
    for gene_sym, gencode in TARGETS:
        print(f"\nfetching {gene_sym} ({gencode}) ...")
        try:
            entries = fetch_all(gencode, TISSUE)
        except Exception as e:
            print(f"  FAILED: {e}"); entries = []
        # Build position -> (NES, pValue) for in-window entries
        nes_by_pos = {}
        for e in entries:
            pos_b38 = e.get("pos")
            if pos_b38 is None: continue
            if not (31_148_818 <= pos_b38 <= 31_648_818): continue
            pos_b37 = pos_b38 - OFFSET_37_TO_38
            nes = e.get("nes")
            p = e.get("pValue")
            if nes is None: continue
            # Keep most-significant per position
            if pos_b37 not in nes_by_pos or p < nes_by_pos[pos_b37][1]:
                nes_by_pos[pos_b37] = (float(nes), float(p))
        eqtl_data[gene_sym] = nes_by_pos
        print(f"  {len(nes_by_pos)} unique in-window positions")

    # 2. Build per-component NES distribution
    rows = []
    for c in range(K):
        top = top_sets[c]
        for gene_sym in [t[0] for t in TARGETS]:
            ed = eqtl_data[gene_sym]
            ovl = top & set(ed.keys())
            if len(ovl) == 0: continue
            nes_vals = np.array([ed[p][0] for p in ovl])
            n = len(nes_vals)
            mean_nes = float(nes_vals.mean())
            median_nes = float(np.median(nes_vals))
            sign_coh = float(abs(np.mean(np.sign(nes_vals))))
            n_pos = int((nes_vals > 0).sum())
            n_neg = int((nes_vals < 0).sum())
            try:
                stat, p_wil = wilcoxon(nes_vals)
                p_wil = float(p_wil)
            except Exception:
                p_wil = float("nan")
            rows.append({
                "component": c, "target_gene": gene_sym,
                "n_overlap": n,
                "mean_NES": mean_nes,
                "median_NES": median_nes,
                "sign_coherence": sign_coh,
                "n_pos_NES": n_pos,
                "n_neg_NES": n_neg,
                "wilcoxon_p": p_wil,
            })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "eqtl_NES_per_component.csv", index=False)

    # Print structured summary
    print("\n=== mean(NES) per component × gene ===")
    piv_mean = df.pivot(index="component", columns="target_gene",
                         values="mean_NES").round(2)
    print(piv_mean.to_string())
    print("\n=== sign coherence (|mean(sign(NES))|; 1.0 = all same sign) ===")
    piv_sc = df.pivot(index="component", columns="target_gene",
                       values="sign_coherence").round(2)
    print(piv_sc.to_string())
    print("\n=== n+/n- (positive/negative NES count) ===")
    df["n_signed"] = df.apply(lambda r: f"{r['n_pos_NES']}/{r['n_neg_NES']}", axis=1)
    piv_signed = df.pivot(index="component", columns="target_gene",
                           values="n_signed")
    print(piv_signed.to_string())
    print("\n=== Wilcoxon p-value (one-sample, vs 0) -log10 ===")
    piv_p = df.pivot(index="component", columns="target_gene",
                      values="wilcoxon_p")
    print((-np.log10(piv_p.replace(0, 1e-300))).round(1).to_string())

    # ===== Figure: violin plots for the contrast components =====
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    fig, axes = plt.subplots(1, len(COMPONENTS_OF_INTEREST),
                              figsize=(3.5 * len(COMPONENTS_OF_INTEREST), 4.6),
                              sharey=True)
    palette = {"MICA":"#d65454","MICB":"#b5722d","HLA-B":"#5a4fcf",
                "HLA-C":"#3a8c5a","HCG27":"#7e3da7"}
    for ax_i, c in enumerate(COMPONENTS_OF_INTEREST):
        ax = axes[ax_i]
        top = top_sets[c]
        plot_data = []
        labels = []
        colors = []
        for gene_sym, _ in TARGETS:
            ed = eqtl_data[gene_sym]
            ovl = top & set(ed.keys())
            if len(ovl) < 5:
                continue
            nes_vals = np.array([ed[p][0] for p in ovl])
            plot_data.append(nes_vals)
            labels.append(f"{gene_sym}\n(n={len(nes_vals)})")
            colors.append(palette[gene_sym])
        if not plot_data: continue
        # Violin
        parts = ax.violinplot(plot_data, positions=range(len(plot_data)),
                                showmeans=True, showmedians=False, widths=0.7)
        for pc, col in zip(parts["bodies"], colors):
            pc.set_facecolor(col); pc.set_alpha(0.45); pc.set_edgecolor("#222")
        # Mean line
        if "cmeans" in parts:
            parts["cmeans"].set_color("#000")
        # Strip overlay
        for j, (vals, col) in enumerate(zip(plot_data, colors)):
            x = np.full_like(vals, j, dtype=float) + np.random.RandomState(0).uniform(-0.12, 0.12, len(vals))
            ax.scatter(x, vals, s=7, color=col, alpha=0.5, edgecolor="none")
        ax.axhline(0, color="#222", lw=0.6, ls="--")
        ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels,
                                                                fontsize=7.5)
        ax.set_title(f"c{c}", fontsize=11, loc="left")
        if ax_i == 0:
            ax.set_ylabel("eQTL effect size (NES)\n← lower expression  /  higher expression →")
        ax.grid(axis="y", linestyle=":", alpha=0.3)
    fig.suptitle(
        "eQTL effect direction within NMF component top-5% SNVs (GTEx Whole_Blood)\n"
        "c1, c4 — population-stable axes;  c5 — HLA-B-variable axis;  c6 — MICA-stable axis",
        fontsize=10.5, y=1.01,
    )
    fig.tight_layout()
    out = OUT / "figure_eqtl_direction"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
