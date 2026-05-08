"""GTEx Liver replication of c5 (HLA-B-variable axis) eQTL direction.

For each NMF component (k=8):
  fetch GTEx Liver eQTLs for MICA / MICB / HLA-B / HLA-C / HCG27;
  for top-5% SNVs, compute mean NES, sign coherence, n+/n-.

Compare to Whole_Blood result.

Output:
  results/liver_NES_per_component.csv
  results/figure_liver_replication.{pdf,png}
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
EQTL_NES_BLOOD = REPO / "results/eqtl_NES_per_component.csv"
OUT = REPO / "results"

OFFSET_37_TO_38 = 32_223
TISSUE = "Liver"
TARGETS = [
    ("MICA",   "ENSG00000204520.12"),
    ("MICB",   "ENSG00000204516.9"),
    ("HLA-B",  "ENSG00000234745.10"),
    ("HLA-C",  "ENSG00000204525.16"),
    ("HCG27",  "ENSG00000206344.7"),
]
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
    K = h_mat.shape[0]
    top_n = max(1, int(NMF_TOP_FRAC * len(snv_pos)))

    # Top-5% per component
    top_sets = {}
    for c in range(K):
        idx = np.argsort(h_mat[c])[-top_n:]
        top_sets[c] = set(int(p) for p in snv_pos[idx])

    # Fetch Liver eQTLs per gene
    eqtl_data = {}
    for gene_sym, gencode in TARGETS:
        print(f"\nfetching {gene_sym} ({gencode}) Liver eQTLs ...")
        try:
            entries = fetch_all(gencode, TISSUE)
        except Exception as e:
            print(f"  FAILED: {e}"); entries = []
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
        eqtl_data[gene_sym] = nes_by_pos
        print(f"  {len(nes_by_pos)} unique in-window positions")

    # Per-component analysis
    rows = []
    for c in range(K):
        top = top_sets[c]
        for gene_sym in [t[0] for t in TARGETS]:
            ed = eqtl_data[gene_sym]
            ovl = top & set(ed.keys())
            if len(ovl) == 0:
                rows.append({
                    "component": c, "tissue": TISSUE, "gene": gene_sym,
                    "n_eqtl_in_window": len(ed),
                    "n_overlap": 0,
                    "mean_NES": float("nan"), "median_NES": float("nan"),
                    "n_pos": 0, "n_neg": 0, "sign_coherence": float("nan"),
                })
                continue
            nes_arr = np.array([ed[p][0] for p in ovl])
            n_pos = int((nes_arr > 0).sum()); n_neg = int((nes_arr < 0).sum())
            rows.append({
                "component": c, "tissue": TISSUE, "gene": gene_sym,
                "n_eqtl_in_window": len(ed),
                "n_overlap": len(ovl),
                "mean_NES": float(nes_arr.mean()),
                "median_NES": float(np.median(nes_arr)),
                "n_pos": n_pos, "n_neg": n_neg,
                "sign_coherence": float(abs(np.mean(np.sign(nes_arr)))),
            })
    df_liver = pd.DataFrame(rows)
    df_liver.to_csv(OUT / "liver_NES_per_component.csv", index=False)
    print(f"\nwrote {OUT/'liver_NES_per_component.csv'}")

    # Compare with Whole_Blood
    df_blood = pd.read_csv(EQTL_NES_BLOOD)
    df_blood = df_blood.rename(columns={"target_gene":"gene"})
    df_blood["tissue"] = "Whole_Blood"
    keep = ["component","tissue","gene","n_overlap","mean_NES","sign_coherence","n_pos_NES","n_neg_NES"]
    df_blood_keep = df_blood[["component","tissue","gene","n_overlap","mean_NES",
                                "sign_coherence","n_pos_NES","n_neg_NES"]].rename(
        columns={"n_pos_NES":"n_pos","n_neg_NES":"n_neg"})
    df_combined = pd.concat([
        df_blood_keep,
        df_liver[["component","tissue","gene","n_overlap","mean_NES",
                   "sign_coherence","n_pos","n_neg"]],
    ], ignore_index=True, sort=False)

    # Per-component, per-gene compare
    print("\n=== Whole_Blood vs Liver: mean NES (c5 = HLA-B-variable axis) ===")
    print(df_combined[df_combined["component"] == 5].to_string(index=False))

    print("\n=== mean NES per (component × gene × tissue) — pivot ===")
    piv = df_combined.pivot_table(index=["gene","tissue"], columns="component",
                                     values="mean_NES")
    print(piv.round(2).to_string())

    # ===== Figure: Whole_Blood vs Liver mean NES per component (c5 highlight) =====
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharex=False)
    palette = {"MICA":"#d65454","MICB":"#b5722d","HLA-B":"#5a4fcf",
                "HLA-C":"#3a8c5a","HCG27":"#7e3da7"}

    for c, ax in zip(range(K), axes.flat):
        sub_b = df_combined[(df_combined["component"] == c) &
                              (df_combined["tissue"] == "Whole_Blood")]
        sub_l = df_combined[(df_combined["component"] == c) &
                              (df_combined["tissue"] == "Liver")]
        # Bar plot per gene
        genes = [g for g, _ in TARGETS]
        x = np.arange(len(genes))
        for offset, sub, hatch, lab in [(-0.18, sub_b, None, "Whole Blood"),
                                          (+0.18, sub_l, "//", "Liver")]:
            yv = []; ya = []
            for g in genes:
                row = sub[sub["gene"] == g]
                if len(row) and not np.isnan(row.iloc[0]["mean_NES"]):
                    yv.append(row.iloc[0]["mean_NES"])
                    ya.append(row.iloc[0]["n_overlap"])
                else:
                    yv.append(np.nan); ya.append(0)
            ax.bar(x + offset, yv, width=0.34,
                    color=[palette[g] for g in genes], hatch=hatch,
                    edgecolor="#222", linewidth=0.5, alpha=0.85, label=lab)
            for i, (v, n) in enumerate(zip(yv, ya)):
                if not np.isnan(v):
                    ax.text(x[i] + offset, v + (0.02 if v >= 0 else -0.02),
                             f"n={n}", ha="center",
                             va="bottom" if v >= 0 else "top",
                             fontsize=6.2, color="#444")
        ax.axhline(0, color="#222", lw=0.5)
        ax.set_xticks(x); ax.set_xticklabels(genes, rotation=30, ha="right",
                                                fontsize=7.5)
        ax.set_title(f"c{c}", fontsize=10, loc="left")
        ax.set_ylabel("mean NES")
        ax.grid(axis="y", linestyle=":", alpha=0.3)
        if c == 0:
            ax.legend(fontsize=7.5, loc="best")
    fig.suptitle(
        "Whole_Blood vs Liver: c5 (HLA-B-variable) direction replication\n"
        "(c5 = HLA-B↑ HLA-C↓ MICB↑ in blood — does Liver agree?)",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    out = OUT / "figure_liver_replication"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
