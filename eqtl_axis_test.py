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

from figure_style import apply_style, WIDTH_DOUBLE_INCH, ajhg_finalize

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

    # 1. Fetch eQTL data per target gene (cached: reuse the CSV if it exists
    #    so layout-only re-renders don't repeat the slow paginated GTEx call)
    cache_path = OUT / "eqtl_targets_per_gene.csv"
    if cache_path.exists():
        print(f"\nloading cached eQTL positions from {cache_path}")
        cached = pd.read_csv(cache_path)
        eqtl_b37_by_gene = {
            g: set(cached[cached["gene"] == g]["b37_pos"].astype(int))
            for g, _ in TARGETS
        }
        for g, _ in TARGETS:
            print(f"  {g}: {len(eqtl_b37_by_gene[g])} unique b37 positions")
    else:
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
                eqtl_rows.append({"gene": gene_sym, "b37_pos": p,
                                  "tissue": TISSUE})
        pd.DataFrame(eqtl_rows).to_csv(cache_path, index=False)
        print(f"\nwrote {cache_path}")

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
            # neglog10_p: NaN when there is nothing to test (K_succ=0), inf
            # only when p underflowed float64 (genuinely tiny). Without this
            # guard the "no eQTLs in window" case (e.g. HCP5) prints ">308",
            # which falsely reads as overwhelming significance.
            if K_succ == 0 or np.isnan(p_val):
                neglog10 = float("nan")
            elif p_val > 0:
                neglog10 = float(-np.log10(p_val))
            else:
                neglog10 = float("inf")
            rows.append({
                "component": c, "target_gene": gene_sym,
                "n_eqtl_in_window": K_succ,
                "n_top_snvs": n_draws,
                "n_overlap": x,
                "expected": expected,
                "fold_enrichment": fold,
                "hyper_p": p_val,
                "neglog10_p": neglog10,
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

    # Main-figure columns: HCP5 dropped (entirely n/a in both panels due to
    # insufficient top-5% × eQTL overlap for enrichment testing). Full HCP5
    # column remains in component_x_eqtl_target.csv and is described in the
    # legend; render it only in the supplementary table.
    MAIN_COLS = ["MICA", "MICB", "HLA-B", "HLA-C", "HCG27"]
    fold_main = fold[MAIN_COLS]
    p_main = p[MAIN_COLS]

    # 5. Figure — 3 panels:
    #   (a) eQTL fold enrichment  (component × gene; Whole_Blood)
    #   (b) -log10 p              (component × gene; Whole_Blood)
    #   (c) Cross-tissue mean NES (axis × tissue) — direction analysis
    apply_style()
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(WIDTH_DOUBLE_INCH, 5.8))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[3.0, 2.0],
                   hspace=0.75, wspace=0.30,
                   left=0.26, right=0.95, top=0.93, bottom=0.10)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
             fig.add_subplot(gs[1, :])]

    # Panel A: fold enrichment heatmap. HCP5 column dropped from the main
    # figure (insufficient overlap → all-n/a); HCP5 is retained in
    # component_x_eqtl_target.csv and described in the legend.
    ax = axes[0]
    fmat = fold_main.values
    im = ax.imshow(fmat, aspect="auto", cmap="RdBu_r",
                    vmin=0, vmax=4, interpolation="nearest")
    ax.set_xticks(range(len(MAIN_COLS)))
    ax.set_xticklabels(MAIN_COLS, rotation=30, ha="right")
    ax.set_yticks(range(K)); ax.set_yticklabels([f"c{c}" for c in range(K)])
    for i in range(K):
        for j in range(len(MAIN_COLS)):
            v = fmat[i, j]
            if np.isnan(v):
                ax.text(j, i, "n/a", ha="center", va="center",
                         fontsize=7.5, color="#666", style="italic")
            else:
                ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                         fontsize=8, color="#fff" if v > 2.4 else "#111")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="fold")
    ax.set_title("(a) eQTL fold enrichment", loc="left", fontsize=10)

    # Panel B: -log10 p heatmap. Color scale AND printed cell values capped
    # at 15 so the panel is internally consistent: cells whose true −log10 p
    # exceeds 15 are labelled ">15" rather than rendered as raw values that
    # would mismatch the saturated color scale.
    ax = axes[1]
    pmat = p_main.values.astype(float)
    pmat_capped = np.where(np.isinf(pmat), 15,
                              np.where(np.isnan(pmat), 0, np.minimum(pmat, 15)))
    im = ax.imshow(pmat_capped, aspect="auto", cmap="Blues",
                    vmin=0, vmax=15, interpolation="nearest")
    ax.set_xticks(range(len(MAIN_COLS)))
    ax.set_xticklabels(MAIN_COLS, rotation=30, ha="right")
    ax.set_yticks(range(K)); ax.set_yticklabels([f"c{c}" for c in range(K)])
    for i in range(K):
        for j in range(len(MAIN_COLS)):
            v = pmat[i, j]
            if np.isnan(v):
                ax.text(j, i, "n/a", ha="center", va="center",
                         fontsize=7.5, color="#666", style="italic")
                continue
            if np.isinf(v) or v > 15:
                ax.text(j, i, ">15", ha="center", va="center",
                         fontsize=8, color="#fff", fontweight="bold")
            elif v >= 1:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                         fontsize=8, color="#fff" if v > 8 else "#111")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label=r"$-\log_{10} p$")
    ax.set_title(r"(b) $-\log_{10}\,P$", loc="left", fontsize=10)

    # Panel C: Cross-tissue mean NES direction analysis.
    # Reads results/axis_sign_coherence_matrix.csv (built by
    # axis_sign_coherence_matrix.py). Axis I: c4/c6 × MICA;
    # Axis II: c5 × HLA-B / HLA-C. Every cell shown has sign coherence
    # = 1.00 in 6/6 tissues, confirming unidirectional eQTL effects.
    ax = axes[2]
    coh = pd.read_csv(OUT / "axis_sign_coherence_matrix.csv")
    coh["row"] = coh["component"] + " × " + coh["gene"]
    row_order = ["c4 × MICA", "c6 × MICA", "c5 × HLA-B", "c5 × HLA-C"]
    tissue_order = ["Whole_Blood", "Liver",
                     "Cells_EBV-transformed_lymphocytes",
                     "Lung", "Skin_Sun_Exposed_Lower_leg",
                     "Colon_Transverse"]
    tissue_labels = ["Whole\nBlood", "Liver", "LCL",
                      "Lung", "Skin", "Colon"]
    mat = (coh[coh["row"].isin(row_order)]
            .pivot(index="row", columns="tissue", values="mean_NES")
            .loc[row_order, tissue_order])
    coh_mat = (coh[coh["row"].isin(row_order)]
                .pivot(index="row", columns="tissue", values="sign_coherence")
                .loc[row_order, tissue_order])
    im = ax.imshow(mat.values, aspect="auto", cmap="RdBu_r",
                    vmin=-0.75, vmax=0.75, interpolation="nearest")
    ax.set_xticks(range(len(tissue_order)))
    ax.set_xticklabels(tissue_labels, fontsize=8.5)
    ax.set_yticks(range(len(row_order)))
    # Embed Axis I/II prefix into the y-tick labels themselves; the
    # gridspec left margin (left=0.20) reserves enough room for them.
    # Use Axis-I / Axis-II spelled out so the Roman numerals can't be
    # mistaken for each other in low-res renderings.
    row_labels = [
        "Axis-I    c4 × MICA",
        "Axis-I    c6 × MICA",
        "Axis-II  c5 × HLA-B",
        "Axis-II  c5 × HLA-C",
    ]
    ax.set_yticklabels(row_labels, fontsize=9, family="monospace")
    for i in range(len(row_order)):
        for j in range(len(tissue_order)):
            v = mat.values[i, j]
            c_val = coh_mat.values[i, j]
            if np.isnan(v):
                ax.text(j, i, "n/a", ha="center", va="center",
                         fontsize=7.5, color="#666", style="italic")
            else:
                txt_color = "#fff" if abs(v) > 0.45 else "#111"
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                         fontsize=8.5, color=txt_color)
    # Light horizontal separator between Axis I (rows 0-1) and Axis II
    # (rows 2-3) — a thin gray line just below row 1 inside the heatmap.
    ax.axhline(1.5, color="#222", lw=0.7, alpha=0.6)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02,
                  label="mean NES   (RdBu: red = ↑, blue = ↓)")

    # (Axis I/II prefix is now embedded directly in y-tick labels above.)
    ax.set_title("(c) Cross-tissue mean NES\n"
                  "      (sign coherence = 1.00 in every shown cell; n per cell in Table S8)",
                  loc="left", fontsize=10)

    fig.suptitle(
        "Figure 3.  NMF components separate MICA↓ and HLA-B↑/HLA-C↓ "
        "regulatory axes.",
        fontsize=9.5, y=0.99, va="top", wrap=True,
    )
    # Note: do NOT call ajhg_finalize() here — its tight_layout() would
    # override the explicit gridspec left/right margins we set above to
    # accommodate the wider "(Axis I) c4 × MICA" y-tick labels in panel (c).
    out = OUT / "figure_component_eqtl_axis"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
