"""Test whether MICA-stable axis (c4, c6) and HLA-B-variable axis (c3, c5)
correspond to distinct expression-regulation signals.

Approach:
  1. Fetch GTEx Whole_Blood eQTL SNVs for 6 target genes:
        MICA, MICB, HLA-B, HLA-C, HCP5, HCG27
     Each yields a set of significant eQTL positions (GRCh38).
  2. Convert to GRCh37 via offset (-32,223 for chr6 in this region).
  3. For each NMF k=8 component, compute:
        - top 5% SNV set (signature)
        - overlap with each target-gene eQTL set
        - hypergeometric enrichment p-value, fold
  4. Report per-component eQTL profile; test directional hypothesis.

Outputs:
  results/eqtl_targets_per_gene.csv   (b37 positions per target gene)
  results/component_x_eqtl_target.csv (component × gene enrichment)
  results/figure_component_eqtl_axis.{pdf,png}
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
from scipy.stats import hypergeom

REPO = Path(__file__).resolve().parent
H_PARQUET = REPO / "results/nmf_H_26_k8.parquet"
OUT = REPO / "results"

OFFSET_37_TO_38 = 32_223   # chr6 GRCh37 -> GRCh38 lift in this region

TARGETS = [
    ("MICA",   "ENSG00000204520.12"),
    ("MICB",   "ENSG00000204516.9"),
    ("HLA-B",  "ENSG00000234745.10"),
    ("HLA-C",  "ENSG00000204525.16"),
    ("HCP5",   "ENSG00000206337.10"),
    ("HCG27",  "ENSG00000206344.7"),
]
TISSUE = "Whole_Blood"
NMF_TOP_FRAC = 0.05


def fetch_all(gencode: str, tissue: str):
    """Fetch all paginated eQTLs for (gene, tissue). Returns list of dicts."""
    out = []
    page = 0
    items_per_page = 2000
    while True:
        url = (f"https://gtexportal.org/api/v2/association/singleTissueEqtl"
                f"?gencodeId={gencode}&tissueSiteDetailId={tissue}"
                f"&datasetId=gtex_v8&itemsPerPage={items_per_page}&page={page}")
        req = urllib.request.Request(url, headers={
            "Accept": "application/json", "User-Agent": "Mozilla/5.0",
        })
        with urllib.request.urlopen(req, timeout=60) as r:
            data = json.loads(r.read().decode())
        out.extend(data.get("data", []))
        info = data.get("paging_info", {})
        n_pages = info.get("numberOfPages", 1)
        if page + 1 >= n_pages:
            break
        page += 1
        time.sleep(0.3)
    return out


def main():
    sys.stdout.reconfigure(line_buffering=True)
    H = pd.read_parquet(H_PARQUET)
    snv_pos_b37 = np.array([int(c) for c in H.columns])
    snv_pos_b38 = snv_pos_b37 + OFFSET_37_TO_38
    print(f"NMF: {H.shape[0]} components × {H.shape[1]} SNVs")

    # 1. Fetch eQTL data per target gene
    eqtl_rows = []
    eqtl_b37_by_gene = {}
    for gene_sym, gencode in TARGETS:
        print(f"\nfetching {gene_sym} ({gencode}) eQTLs in {TISSUE} ...")
        t0 = time.time()
        try:
            entries = fetch_all(gencode, TISSUE)
        except Exception as e:
            print(f"  FAILED: {e}")
            entries = []
        print(f"  {len(entries)} entries in {time.time()-t0:.0f}s")
        positions_b38 = sorted({int(e.get("pos", -1)) for e in entries
                                  if e.get("pos") is not None})
        # Filter to our window (b38: 31,148,818 - 31,648,818)
        in_window = [p for p in positions_b38
                     if 31_148_818 <= p <= 31_648_818]
        # Convert to b37 (subtract offset)
        positions_b37 = set(p - OFFSET_37_TO_38 for p in in_window)
        eqtl_b37_by_gene[gene_sym] = positions_b37
        print(f"  in-window (b38): {len(in_window)}; "
               f"unique b37 positions: {len(positions_b37)}")
        for p in positions_b37:
            eqtl_rows.append({"gene": gene_sym, "b37_pos": p, "tissue": TISSUE})
    pd.DataFrame(eqtl_rows).to_csv(OUT / "eqtl_targets_per_gene.csv", index=False)
    print(f"\nwrote {OUT/'eqtl_targets_per_gene.csv'}")

    # 2. NMF top-5% SNV set per component
    h_mat = H.values
    K = h_mat.shape[0]
    top_frac = NMF_TOP_FRAC
    top_sets = {}
    for c in range(K):
        h = h_mat[c]
        top_n = max(1, int(top_frac * len(h)))
        top_idx = np.argsort(h)[-top_n:]
        top_sets[c] = set(int(p) for p in snv_pos_b37[top_idx])

    # 3. Universe = all SNVs in NMF
    universe = set(int(p) for p in snv_pos_b37)
    N_uni = len(universe)

    # 4. Hypergeometric enrichment per (component × target gene)
    rows = []
    for c in range(K):
        top = top_sets[c]
        for gene_sym, _ in TARGETS:
            eqtl_set = eqtl_b37_by_gene[gene_sym] & universe
            x = len(top & eqtl_set)
            K_succ = len(eqtl_set)
            n_draws = len(top)
            if K_succ == 0:
                p_val = float("nan"); fold = float("nan"); expected = 0
            else:
                p_val = float(hypergeom.sf(x - 1, N_uni, K_succ, n_draws))
                expected = n_draws * K_succ / N_uni
                fold = x / expected if expected > 0 else float("nan")
            rows.append({
                "component": c, "target_gene": gene_sym,
                "n_eqtl_in_window": K_succ,
                "n_top_snvs": n_draws,
                "n_overlap": x,
                "expected": expected,
                "fold_enrichment": fold,
                "hyper_p": p_val,
                "neglog10_p": -np.log10(p_val) if (p_val and p_val > 0) else float("inf"),
            })
    enrich = pd.DataFrame(rows)
    enrich.to_csv(OUT / "component_x_eqtl_target.csv", index=False)

    # Print summary
    print("\n=== Component × target gene fold enrichment ===")
    fold = enrich.pivot(index="component", columns="target_gene",
                         values="fold_enrichment").round(2)
    cols = [g for g, _ in TARGETS]
    fold = fold[cols]
    print(fold.to_string())

    print("\n=== Component × target gene -log10(p-value) ===")
    p = enrich.pivot(index="component", columns="target_gene",
                      values="neglog10_p").round(1)
    p = p[cols]
    print(p.to_string())

    # 5. Figure
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))

    # Panel A: fold enrichment heatmap
    ax = axes[0]
    fmat = fold.values
    im = ax.imshow(fmat, aspect="auto", cmap="RdBu_r",
                    vmin=0, vmax=4, interpolation="nearest")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols,
                                                          rotation=30, ha="right")
    ax.set_yticks(range(K)); ax.set_yticklabels([f"c{c}" for c in range(K)])
    for i in range(K):
        for j in range(len(cols)):
            v = fmat[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                         fontsize=8, color="#fff" if v > 2.4 else "#111")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="fold")
    ax.set_title("A. NMF component × target-gene eQTL fold enrichment\n"
                  "(top 5% SNVs vs Whole_Blood eQTLs)", loc="left", fontsize=10)

    # Panel B: -log10 p heatmap
    ax = axes[1]
    pmat = p.values.astype(float)
    pmat_capped = np.where(np.isnan(pmat), 0, np.minimum(pmat, 30))
    im = ax.imshow(pmat_capped, aspect="auto", cmap="Blues",
                    vmin=0, vmax=15, interpolation="nearest")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols,
                                                          rotation=30, ha="right")
    ax.set_yticks(range(K)); ax.set_yticklabels([f"c{c}" for c in range(K)])
    for i in range(K):
        for j in range(len(cols)):
            v = pmat[i, j]
            if not np.isnan(v) and v >= 1:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                         fontsize=8, color="#fff" if v > 8 else "#111")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label=r"$-\log_{10} p$")
    ax.set_title("B. Hypergeometric significance (capped at 30)",
                  loc="left", fontsize=10)

    fig.suptitle(
        "Do NMF components correspond to distinct eQTL signals?\n"
        "(rs2596542-T anchor-carrier branches, GTEx v8 Whole_Blood)",
        fontsize=11, y=1.04,
    )
    fig.tight_layout()
    out = OUT / "figure_component_eqtl_axis"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
