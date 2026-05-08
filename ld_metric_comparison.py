"""(e) Compare flip detection across LD metrics.

For 8 metrics, define "flip" between super-pops:
  Signed metrics (sign change between super-pops):
    signed_r       — Pearson r (already in CSV)
    D              — covariance p11 - p_a*p_b
    D_prime        — Lewontin's D' (D / D_max)
    Delta          — P(A|B) − P(B|A)
    C              — Δ − M (structural component)

  Magnitude metrics (no sign — use IQR / range as proxy for "differs"):
    r2             — r² (unsigned)
    abs_r          — |r| (unsigned)
    M              — p_a − p_b (signed but trivially flips with anchor coding)

For each NMF component (k=8) × metric, hypergeom test of top-5% SNV
overlap with the corresponding flip set (any super-pop pair flip in
sign for signed metrics; any super-pop pair difference > 0.2 for
unsigned proxies — i.e., trying to give unsigned metrics a fair shot).

Output:
  results/ld_metric_flip_per_component.csv
  results/figure_ld_metric_comparison.{pdf,png}
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

REPO = Path(__file__).resolve().parent
H_PARQUET = REPO / "results/nmf_H_26_k8.parquet"
RPP_PARQUET = REPO / "results/region_per_pop_26.parquet"
OUT = REPO / "results"

ANCHOR_POS = 31_366_595
SUPER_POPS = ["AFR", "EUR", "EAS", "SAS", "AMR"]
NMF_TOP_FRAC = 0.05
UNSIGNED_DIFF_THRESHOLD = 0.2   # for unsigned metrics: "flip" = max - min >= 0.2


def derive_per_pop_metrics(df_long: pd.DataFrame) -> pd.DataFrame:
    """Compute D, D_prime, Delta, |r|, r^2, abs_r, signed_r, M, C per (partner_pos, pop)."""
    n11 = df_long["n11"].astype(float)
    p_a = df_long["p_a"].astype(float)
    p_b = df_long["p_b"].astype(float)
    n   = df_long["n"].astype(float)
    p_ab = n11 / n
    q_a = 1 - p_a; q_b = 1 - p_b
    D = p_ab - p_a * p_b
    # D_max: depends on D sign
    D_max = np.where(D >= 0,
                       np.minimum(p_a * q_b, q_a * p_b),
                       np.minimum(p_a * p_b, q_a * q_b))
    D_prime = np.where(D_max > 0, D / D_max, np.nan)
    denom = p_a * q_a * p_b * q_b
    r2 = np.where(denom > 0, (D ** 2) / denom, np.nan)
    abs_r = np.sqrt(np.where(np.isnan(r2), np.nan, np.maximum(r2, 0)))
    signed_r = np.where((D < 0) & ~np.isnan(abs_r), -abs_r, abs_r)
    P_AgB = np.where((n11 + (n - n11 - df_long["n10"] - df_long["n01"]) * 0  # n01 = ?
                        + df_long["n01"]) > 0, np.nan, np.nan)
    # Easier: use existing C / Delta computation pattern
    n10 = df_long["n10"].astype(float)
    n01 = df_long["n01"].astype(float)
    P_AgB = np.where((n11 + n01) > 0, n11 / (n11 + n01), np.nan)
    P_BgA = np.where((n11 + n10) > 0, n11 / (n11 + n10), np.nan)
    Delta = P_AgB - P_BgA
    M = p_a - p_b
    C = Delta - M

    df = df_long.copy()
    df["D"] = D
    df["D_prime"] = D_prime
    df["r2"] = r2
    df["abs_r"] = abs_r
    df["signed_r2"] = signed_r * abs_r  # signed r-squared (rare but useful)
    df["Delta"] = Delta
    df["M_metric"] = M
    df["C_metric"] = C
    return df


def main():
    H = pd.read_parquet(H_PARQUET)
    snv_pos = np.array([int(c) for c in H.columns])
    h_mat = H.values
    K = h_mat.shape[0]
    universe = set(int(p) for p in snv_pos)

    print("loading region_per_pop_26.parquet ...")
    rpp = pd.read_parquet(RPP_PARQUET)
    rpp["super_pop"] = rpp["pop"].apply(
        lambda p: ({s: ssup for ssup in ["AFR","EUR","EAS","SAS","AMR"]
                     for s in [ssup, f"{ssup}_pooled"]}.get(p, None)
                    if not p.endswith("_pooled") else p.split("_")[0]))
    # Actually simpler: just map pop -> super_pop via panel
    pop_to_super = {}
    pmap = {"AFR":["YRI","LWK","GWD","MSL","ESN","ASW","ACB"],
             "EUR":["CEU","FIN","GBR","IBS","TSI"],
             "EAS":["CHB","JPT","CHS","CDX","KHV"],
             "SAS":["GIH","PJL","BEB","STU","ITU"],
             "AMR":["MXL","PUR","CLM","PEL"]}
    for s, ps in pmap.items():
        for p in ps:
            pop_to_super[p] = s
    rpp["super_pop"] = rpp["pop"].map(pop_to_super)
    sub = rpp.dropna(subset=["super_pop"]).copy()

    # Need to compute Delta_user = P(A|B) - P(B|A) (user convention) plus M, C
    # The CSV-format file has columns: n,n11,p_a,p_b,M,C,Delta,r_signed,topology_class
    # Verify: re-check columns
    print("rpp columns:", list(sub.columns))

    # Already have M, C, Delta, r_signed. Compute D, D', |r|, r²
    n11 = sub["n11"].astype(float)
    p_a = sub["p_a"].astype(float)
    p_b = sub["p_b"].astype(float)
    n   = sub["n"].astype(float)
    p_ab = n11 / n
    q_a = 1 - p_a; q_b = 1 - p_b
    D = p_ab - p_a * p_b
    D_max = np.where(D >= 0,
                       np.minimum(p_a * q_b, q_a * p_b),
                       np.minimum(p_a * p_b, q_a * q_b))
    sub["D"] = D
    sub["D_prime"] = np.where(D_max > 0, D / D_max, np.nan)
    sub["r_squared"] = sub["r_signed"].astype(float) ** 2
    sub["abs_r"] = sub["r_signed"].astype(float).abs()

    # Compute median per (partner_pos, super_pop) for each metric
    metrics = ["signed_r", "D", "D_prime", "Delta", "C", "M",
               "r_squared", "abs_r"]
    metric_names = {"signed_r":"r_signed", "D":"D", "D_prime":"D_prime",
                     "Delta":"Delta", "C":"C", "M":"M",
                     "r_squared":"r_squared", "abs_r":"abs_r"}

    # Build per-partner × super-pop median for each metric
    print("computing super-pop medians per metric ...")
    grouped = sub.groupby(["partner_pos","super_pop"])
    medians = grouped[[metric_names[m] for m in metrics]].median().reset_index()

    # Pivot per metric
    flip_sets = {}
    for m in metrics:
        col = metric_names[m]
        pivot = medians.pivot(index="partner_pos", columns="super_pop", values=col)
        if pivot.shape[1] < 5: continue
        mat = pivot[SUPER_POPS].astype(float).values
        is_signed = m in ("signed_r","D","D_prime","Delta","C","M")
        n_flip = np.zeros(len(pivot), dtype=int)
        for i, j in combinations(range(5), 2):
            ci = mat[:, i]; cj = mat[:, j]
            valid = ~(np.isnan(ci) | np.isnan(cj))
            if is_signed:
                # sign flip
                flip = valid & ((np.sign(ci) * np.sign(cj)) < 0)
            else:
                # magnitude difference >= threshold
                flip = valid & (np.abs(ci - cj) >= UNSIGNED_DIFF_THRESHOLD)
            n_flip += flip.astype(int)
        flip_partner = set(int(p) for p, nf in zip(pivot.index, n_flip)
                            if nf >= 1)
        flip_sets[m] = flip_partner & universe
    print(f"flip set sizes (any-pair flip):")
    for m, s in flip_sets.items():
        print(f"  {m}: {len(s)} ({len(s)/len(universe)*100:.1f}%)")

    # Per-component enrichment for each metric
    rows = []
    n_uni = len(universe)
    for c in range(K):
        h = h_mat[c]
        top_n = max(1, int(NMF_TOP_FRAC * h.shape[0]))
        top_idx = np.argsort(h)[-top_n:]
        top_set = set(int(p) for p in snv_pos[top_idx])
        for m in metrics:
            if m not in flip_sets: continue
            ts = flip_sets[m]
            x = len(top_set & ts)
            K_succ = len(ts)
            if K_succ == 0:
                rows.append({"component": c, "metric": m,
                              "fold": float("nan"), "p": float("nan"),
                              "x": x})
                continue
            p_val = float(hypergeom.sf(x - 1, n_uni, K_succ, top_n))
            exp = top_n * K_succ / n_uni
            rows.append({"component": c, "metric": m,
                          "fold": x / exp if exp > 0 else float("nan"),
                          "p": p_val, "x": x, "expected": exp})
    df = pd.DataFrame(rows)
    df["neglog10_p"] = -np.log10(df["p"].replace(0, 1e-300))
    df.to_csv(OUT / "ld_metric_flip_per_component.csv", index=False)

    print("\n=== Fold enrichment per (component × metric) ===")
    fp = df.pivot(index="component", columns="metric", values="fold")
    fp = fp[metrics]
    print(fp.round(2).to_string())
    print("\n=== -log10(p) ===")
    pp = df.pivot(index="component", columns="metric", values="neglog10_p")
    pp = pp[metrics]
    print(pp.round(1).to_string())

    # ===== Figure: heatmap fold enrichment =====
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    fig, ax = plt.subplots(figsize=(11, 4.2))
    fmat = fp.values
    vmax = max(1.6, np.nanmax(fmat) if fmat.size else 1.6)
    im = ax.imshow(fmat, aspect="auto", cmap="RdBu_r",
                    vmin=0.5, vmax=vmax, interpolation="nearest")
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metrics, rotation=20, ha="right")
    ax.set_yticks(range(K)); ax.set_yticklabels([f"c{c}" for c in range(K)])
    for i in range(K):
        for j in range(len(metrics)):
            v = fmat[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                         fontsize=7.5, color="#fff" if abs(v - 1) > 0.4 else "#111")
    fig.colorbar(im, ax=ax, fraction=0.04,
                  label="fold enrichment\n(top-5% SNV in flip set)")
    # Mark signed vs unsigned with x-axis grouping
    ax.axvline(5.5, color="#222", lw=1.2, ls="--", alpha=0.7)
    ax.text(2.5, -1.0, "signed metrics (sign-flip)", ha="center", fontsize=8.5,
             color="#5a4fcf", fontweight="bold")
    ax.text(6.5, -1.0, "unsigned (|Δ|≥0.2)", ha="center", fontsize=8.5,
             color="#888")
    ax.set_title(
        "(e) Flip detection across LD metrics — does the c5 enrichment depend on choice of metric?",
        loc="left", fontsize=10, pad=18,
    )
    fig.tight_layout()
    out = OUT / "figure_ld_metric_comparison"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
