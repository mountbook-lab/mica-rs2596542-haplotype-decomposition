"""External tissue replication of c5 effect direction.

Refetch HLA-B / HLA-C / MICB eQTLs in three additional GTEx tissues:
  - Cells_EBV-transformed_lymphocytes  (LCL — immune-relevant)
  - Lung
  - Skin_Sun_Exposed_Lower_leg

For each tissue × target gene, recompute c5 top-5% NES coherence and
sign distribution. Tests: HLA-B↑, HLA-C↓, MICB↑ direction held in c5
across tissues?

Outputs:
  results/c5_external_tissue_NES.csv
  results/figure_external_tissue_replication.{pdf,png}
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

REPO = Path(__file__).resolve().parent
H_PARQUET = REPO / "results/nmf_H_26_k8.parquet"
OUT = REPO / "results"

OFFSET_37_TO_38 = 32_223
TARGETS = [
    ("HLA-B",  "ENSG00000234745.10"),
    ("HLA-C",  "ENSG00000204525.16"),
    ("MICB",   "ENSG00000204516.9"),
]
TISSUES = [
    "Whole_Blood",                          # baseline (already done)
    "Cells_EBV-transformed_lymphocytes",    # LCL — immune
    "Lung",                                 # immune-active organ
    "Skin_Sun_Exposed_Lower_leg",           # different tissue, sanity check
]
TARGET_COMPONENT = 5
NMF_TOP_FRAC = 0.05


def fetch_all(gencode: str, tissue: str):
    out = []; page = 0
    while True:
        url = (f"https://gtexportal.org/api/v2/association/singleTissueEqtl"
                f"?gencodeId={gencode}&tissueSiteDetailId={tissue}"
                f"&datasetId=gtex_v8&itemsPerPage=2000&page={page}")
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
    snv_pos = np.array([int(c) for c in H.columns])
    h_mat = H.values
    top_n = max(1, int(NMF_TOP_FRAC * len(snv_pos)))
    c5_top_idx = np.argsort(h_mat[TARGET_COMPONENT])[-top_n:]
    c5_top_pos = set(int(p) for p in snv_pos[c5_top_idx])
    print(f"c5 top-5% size: {len(c5_top_pos)}")

    rows = []
    for tissue in TISSUES:
        print(f"\n=== {tissue} ===")
        for gene_sym, gencode in TARGETS:
            print(f"  fetching {gene_sym} ...")
            try:
                entries = fetch_all(gencode, tissue)
            except Exception as e:
                print(f"    FAILED: {e}"); continue
            nes_by_pos = {}
            for e in entries:
                pos_b38 = e.get("pos")
                if pos_b38 is None: continue
                if not (31_148_818 <= pos_b38 <= 31_648_818): continue
                pos_b37 = pos_b38 - OFFSET_37_TO_38
                nes = e.get("nes"); p = e.get("pValue")
                if nes is None: continue
                if pos_b37 not in nes_by_pos or p < nes_by_pos[pos_b37][1]:
                    nes_by_pos[pos_b37] = (float(nes), float(p))
            print(f"    {len(nes_by_pos)} unique in-window positions")
            ovl = c5_top_pos & set(nes_by_pos.keys())
            if not ovl:
                rows.append({
                    "tissue": tissue, "gene": gene_sym,
                    "n_eqtl_in_window": len(nes_by_pos),
                    "n_overlap_c5": 0,
                    "mean_NES": float("nan"), "median_NES": float("nan"),
                    "n_pos": 0, "n_neg": 0, "sign_coherence": float("nan"),
                })
                continue
            nes_arr = np.array([nes_by_pos[p][0] for p in ovl])
            n_pos = int((nes_arr > 0).sum()); n_neg = int((nes_arr < 0).sum())
            rows.append({
                "tissue": tissue, "gene": gene_sym,
                "n_eqtl_in_window": len(nes_by_pos),
                "n_overlap_c5": len(ovl),
                "mean_NES": float(nes_arr.mean()),
                "median_NES": float(np.median(nes_arr)),
                "n_pos": n_pos, "n_neg": n_neg,
                "sign_coherence": float(abs(np.mean(np.sign(nes_arr)))),
            })
            print(f"    overlap={len(ovl)}, mean NES={nes_arr.mean():+.3f}, "
                   f"n+/n-={n_pos}/{n_neg}, coherence={rows[-1]['sign_coherence']:.2f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "c5_external_tissue_NES.csv", index=False)

    # Print pivot: mean NES per tissue × gene
    print("\n=== mean NES of c5 top-5% × tissue × gene ===")
    piv = df.pivot(index="tissue", columns="gene", values="mean_NES")
    piv = piv.reindex(TISSUES)
    print(piv.round(3).to_string())

    print("\n=== sign coherence ===")
    piv_sc = df.pivot(index="tissue", columns="gene", values="sign_coherence")
    piv_sc = piv_sc.reindex(TISSUES)
    print(piv_sc.round(2).to_string())

    print("\n=== n+/n- ===")
    df["n_signed"] = df.apply(lambda r: f"{r['n_pos']}/{r['n_neg']}", axis=1)
    piv_n = df.pivot(index="tissue", columns="gene", values="n_signed")
    piv_n = piv_n.reindex(TISSUES)
    print(piv_n.to_string())

    # ===== Figure: tissue × gene mean NES heatmap with sign coherence overlay =====
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    fig, ax = plt.subplots(figsize=(7, 4))
    nes_mat = piv.values
    sc_mat = piv_sc.values
    n_signed = piv_n.values

    im = ax.imshow(nes_mat, aspect="auto", cmap="RdBu_r",
                    vmin=-0.4, vmax=0.4, interpolation="nearest")
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels([t.replace("_", "\n").replace("-", "-\n") for t in piv.index],
                        fontsize=7.5)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = nes_mat[i, j]; sc = sc_mat[i, j]
            ns = n_signed[i, j] if not pd.isna(n_signed[i, j]) else "0/0"
            if not np.isnan(v):
                txt = f"{v:+.2f}\n{ns}\ncoh={sc:.2f}"
                ax.text(j, i, txt, ha="center", va="center",
                         fontsize=6.8, color="#fff" if abs(v) > 0.25 else "#111")
    fig.colorbar(im, ax=ax, fraction=0.04,
                  label="mean NES (+ ↑ expression)")
    ax.set_title(
        "c5 top-5% NES across tissues — does HLA-B↑ HLA-C↓ MICB↑ replicate?",
        loc="left", fontsize=10,
    )
    fig.tight_layout()
    out = OUT / "figure_external_tissue_replication"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
