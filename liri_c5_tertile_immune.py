"""LIRI-JP c5 tertile boxplot + immune-adjusted regression.

(b) Box plot: c5 tertile (low/mid/high) × log2(HLA-C TPM+1)
    — etiology color, Jonckheere-Terpstra trend test, pairwise tests

(d) Multivariate with immune infiltration proxies (no histological
    purity in Fujimoto S-table1, so use expression-based proxies):
      PTPRC (CD45)  — pan-leukocyte
      CD3D, CD8A    — T cell
      IFNG          — IFN signaling
      B2M, TAP1     — class-I antigen-processing machinery
      HLA-A         — closest co-regulated MHC class I gene
                      (strongest collinearity / strictest control)

To avoid overfitting at n=122, fit four progressively stricter models:
  M1:  c5 + T_stage + viral_status + age + sex      (baseline from prev step)
  M2:  + PTPRC log2(TPM+1)                            (general immune)
  M3:  + immune_PC1   (1st PC of PTPRC/CD3D/CD8A/IFNG/B2M/TAP1)
  M4:  + log2(HLA-A)  (strictest: shared-regulation control)

Output:
  results/liri_c5_tertile_immune.csv
  results/figure_liri_c5_tertile.{pdf,png}
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import (
    kruskal, mannwhitneyu, spearmanr,
)

REPO = Path(__file__).resolve().parent
LIRI_C5 = REPO / "results/liri_c5_score.csv"
CLIN = "/mnt/e/LICA-CN-analysis/LIRI-JP_analysis/data/clinical_data/EGA_clinical_matched.csv"
EXPR = "/mnt/e/LICA-CN-analysis/LIRI-JP_analysis/results/expression_matrix/tpm_gene_named.csv"
OUT = REPO / "results"

IMMUNE_GENES = ["PTPRC", "CD3D", "CD8A", "IFNG", "B2M", "TAP1"]


def jonckheere(values_by_group):
    """Jonckheere-Terpstra one-sided trend statistic + permutation p.
    values_by_group: list of arrays in ordered groups (low → high).
    """
    k = len(values_by_group)
    J = 0
    for i in range(k):
        for j in range(i + 1, k):
            xi = values_by_group[i]; xj = values_by_group[j]
            J += float(((xj[:, None] > xi[None, :]).sum() +
                       0.5 * (xj[:, None] == xi[None, :]).sum()))
    # Two-sided permutation
    rng = np.random.default_rng(0)
    all_v = np.concatenate(values_by_group)
    sizes = [len(v) for v in values_by_group]
    n_perm = 5000
    null = np.zeros(n_perm)
    for p in range(n_perm):
        perm = rng.permutation(all_v)
        groups = []
        s = 0
        for sz in sizes:
            groups.append(perm[s:s+sz]); s += sz
        Jp = 0
        for i in range(k):
            for j in range(i + 1, k):
                xi = groups[i]; xj = groups[j]
                Jp += float(((xj[:, None] > xi[None, :]).sum() +
                            0.5 * (xj[:, None] == xi[None, :]).sum()))
        null[p] = Jp
    p_two = float(2 * min(((null >= J).sum() + 1) / (n_perm + 1),
                          ((null <= J).sum() + 1) / (n_perm + 1)))
    return J, p_two


def main():
    sys.stdout.reconfigure(line_buffering=True)
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(LIRI_C5)
    clin = pd.read_csv(CLIN).rename(columns={"clinical_sample_id": "rk_id"})
    df = df.merge(clin[["rk_id", "age", "gender", "viral_status",
                         "T_stage"]], on="rk_id")
    df = df[df["T_stage"].notna() & df["age"].notna()
             & df["viral_status"].isin(["HBV", "HCV", "NBNC"])].copy()
    print(f"effective n: {len(df)}")

    # Pull immune marker TPMs from EGA matrix
    print(f"loading immune marker expression ...")
    expr = pd.read_csv(EXPR)
    gene_col = expr.columns[0]
    expr_samples = list(expr.columns[1:])
    feats = expr[gene_col].astype(str)

    for g in IMMUNE_GENES:
        m = (feats == g)
        if not m.any():
            print(f"  WARN: {g} not in expression matrix"); continue
        row = expr.loc[m].iloc[0]
        df[g] = df["rk_id"].map(lambda rk: row[rk] if rk in expr_samples else np.nan)
        df[f"log2_{g}"] = np.log2(df[g].astype(float) + 1)

    df["log2_HLA_C"] = np.log2(df["HLA-C"].astype(float) + 1)
    df["log2_HLA_A"] = np.log2(df["HLA-A"].astype(float) + 1)

    # =========================================================
    # (b) Tertile boxplot
    # =========================================================
    cs = df["c5_score"].astype(float)
    df["c5_score_z"] = (cs - cs.mean()) / cs.std(ddof=0)
    q33, q66 = np.quantile(cs, [1/3, 2/3])
    df["c5_tertile"] = pd.cut(cs, bins=[-np.inf, q33, q66, np.inf],
                                labels=["low", "mid", "high"])

    print("\n=== c5 tertile counts ===")
    print(df["c5_tertile"].value_counts().sort_index())

    print("\n=== HLA-C log2 TPM by c5 tertile ===")
    by_t = []
    for t in ["low", "mid", "high"]:
        v = df.loc[df["c5_tertile"] == t, "log2_HLA_C"].dropna().values
        by_t.append(v)
        print(f"  {t}: n={len(v):3d}  median={np.median(v):.3f}  "
              f"mean={np.mean(v):.3f}  std={np.std(v):.3f}")
    kw_stat, kw_p = kruskal(*by_t)
    print(f"  Kruskal-Wallis: stat={kw_stat:.2f}, p={kw_p:.4f}")

    # Jonckheere trend test (5000 permutations, ordered low→high)
    J, p_jt = jonckheere(by_t)
    print(f"  Jonckheere-Terpstra trend (perm): J={J:.0f}, p={p_jt:.4f}")
    # Pairwise low vs high
    u, p_lh = mannwhitneyu(by_t[0], by_t[2], alternative="two-sided")
    print(f"  low vs high MannWhitney: p={p_lh:.4f}")

    # Spearman rho on continuous c5 (more powerful than tertile)
    rs, ps = spearmanr(cs, df["log2_HLA_C"])
    print(f"  Spearman ρ (continuous c5 vs log2 HLA-C): {rs:+.3f}  p={ps:.4f}")

    # =========================================================
    # (d) Immune-adjusted regression
    # =========================================================
    # Build immune_PC1 from 6 markers (z-score → PCA)
    Z = df[[f"log2_{g}" for g in IMMUNE_GENES]].copy()
    for c in Z.columns:
        Z[c] = (Z[c] - Z[c].mean()) / Z[c].std(ddof=0)
    # SVD on rows
    Zv = Z.dropna().values
    U, S, Vt = np.linalg.svd(Zv, full_matrices=False)
    pc_loadings = Vt[0]
    print(f"\n  immune PC1 explained var: {S[0]**2 / (S**2).sum() * 100:.1f}%")
    print(f"  loadings: {dict(zip(IMMUNE_GENES, pc_loadings.round(3).tolist()))}")
    df["immune_PC1"] = (Z.values @ pc_loadings)

    # Encode covariates
    df["sex_M"] = (df["gender"] == "M").astype(int)
    df["viral_HCV"]  = (df["viral_status"] == "HCV").astype(int)
    df["viral_NBNC"] = (df["viral_status"] == "NBNC").astype(int)

    Y = df["log2_HLA_C"]

    models = {
        "M1: + T + viral + age + sex": ["c5_score_z", "T_stage", "viral_HCV",
                                          "viral_NBNC", "age", "sex_M"],
        "M2: + PTPRC":                  ["c5_score_z", "T_stage", "viral_HCV",
                                          "viral_NBNC", "age", "sex_M",
                                          "log2_PTPRC"],
        "M3: + immune_PC1":             ["c5_score_z", "T_stage", "viral_HCV",
                                          "viral_NBNC", "age", "sex_M",
                                          "immune_PC1"],
        "M4: + log2(HLA-A) [strict]":   ["c5_score_z", "T_stage", "viral_HCV",
                                          "viral_NBNC", "age", "sex_M",
                                          "log2_HLA_A"],
    }
    rows = []
    for label, covs in models.items():
        sub = df[[Y.name] + covs].dropna()
        X1 = sm.add_constant(sub[covs])
        m = sm.OLS(sub[Y.name], X1).fit()
        rows.append({
            "model": label, "n": int(m.nobs),
            "c5_beta": float(m.params["c5_score_z"]),
            "c5_se":   float(m.bse["c5_score_z"]),
            "c5_p":    float(m.pvalues["c5_score_z"]),
            "extra_beta": float(m.params[covs[-1]]) if label != "M1: + T + viral + age + sex" else float("nan"),
            "extra_p":    float(m.pvalues[covs[-1]]) if label != "M1: + T + viral + age + sex" else float("nan"),
            "extra_var":  covs[-1] if label != "M1: + T + viral + age + sex" else "",
            "r2":         float(m.rsquared),
        })
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "liri_c5_tertile_immune.csv", index=False)

    print("\n=== Immune-adjusted regression: c5 → log2(HLA-C) ===")
    print(res.to_string(index=False, formatters={
        "c5_beta": "{:+.3f}".format, "c5_se": "{:.3f}".format,
        "c5_p":    "{:.4f}".format,
        "extra_beta": "{:+.3f}".format, "extra_p": "{:.4f}".format,
        "r2":      "{:.3f}".format,
    }))

    # =========================================================
    # Figure
    # =========================================================
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    fig = plt.figure(figsize=(12, 5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.2])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    pal_v = {"HBV":"#d65454", "HCV":"#3a8c5a", "NBNC":"#5a4fcf"}
    rng = np.random.default_rng(0)

    # Box (no fill) + jittered scatter colored by viral_status
    box_data = by_t
    bp = ax1.boxplot(box_data, positions=[0, 1, 2], widths=0.5, patch_artist=False,
                       showfliers=False, medianprops=dict(color="#222", lw=1.5),
                       whiskerprops=dict(color="#888", lw=0.7),
                       capprops=dict(color="#888", lw=0.7),
                       boxprops=dict(color="#888", lw=0.7))
    for i, t in enumerate(["low", "mid", "high"]):
        sub = df[df["c5_tertile"] == t].dropna(subset=["log2_HLA_C", "viral_status"])
        for vs in ["HBV", "HCV", "NBNC"]:
            ssub = sub[sub["viral_status"] == vs]
            if len(ssub) == 0: continue
            xj = i + rng.uniform(-0.15, 0.15, len(ssub))
            ax1.scatter(xj, ssub["log2_HLA_C"], s=22, color=pal_v[vs],
                         alpha=0.7, edgecolor="#222", linewidth=0.4,
                         label=vs if i == 0 else None)
    ax1.set_xticks([0, 1, 2])
    ax1.set_xticklabels([f"low\nn={len(by_t[0])}",
                          f"mid\nn={len(by_t[1])}",
                          f"high\nn={len(by_t[2])}"])
    ax1.set_xlabel("c5 score tertile")
    ax1.set_ylabel(r"$\log_2$(HLA-C TPM + 1)")
    ax1.text(0.02, 0.98,
              f"Jonckheere-Terpstra trend p = {p_jt:.3f}\n"
              f"low vs high MWU p = {p_lh:.3f}\n"
              f"Spearman ρ (continuous c5) = {rs:+.2f}, p = {ps:.3f}",
              transform=ax1.transAxes, fontsize=8.5, va="top",
              bbox=dict(boxstyle="round,pad=0.35", fc="#f8f8f8",
                         ec="#aaa", lw=0.5))
    ax1.legend(title="viral status", fontsize=7.5, title_fontsize=8,
                 loc="lower left", frameon=False)
    ax1.set_title("(b) HLA-C expression by c5 score tertile",
                    loc="left", fontsize=10.5)
    ax1.grid(axis="y", linestyle=":", alpha=0.3)

    # Forest plot of c5 coefficient across models
    labels = list(reversed([r["model"] for _, r in res.iterrows()]))
    betas  = list(reversed([r["c5_beta"] for _, r in res.iterrows()]))
    ses    = list(reversed([r["c5_se"]   for _, r in res.iterrows()]))
    ps     = list(reversed([r["c5_p"]    for _, r in res.iterrows()]))
    ys = np.arange(len(labels))
    for y, b, s, p in zip(ys, betas, ses, ps):
        lo, hi = b - 1.96*s, b + 1.96*s
        col = "#a33" if p < 0.05 else "#888"
        ax2.errorbar([b], [y], xerr=[[b - lo], [hi - b]], fmt="o",
                       color=col, ecolor=col, capsize=4, ms=8, lw=1.5,
                       markeredgecolor="#222", markeredgewidth=0.6)
        ax2.text(hi + 0.02, y,
                  f"  β={b:+.3f}  p={p:.3f}",
                  va="center", fontsize=8.5,
                  color="#222" if p < 0.05 else "#444")
    ax2.axvline(0, color="#222", lw=0.6, ls="-")
    ax2.set_yticks(ys); ax2.set_yticklabels(labels, fontsize=8.5)
    ax2.set_xlabel(r"c5_score (z) $\beta$ on $\log_2$(HLA-C TPM+1) [95% CI]")
    ax2.set_title("(d) c5 → HLA-C robust to immune & class-I co-regulation",
                    loc="left", fontsize=10.5)
    ax2.set_xlim(-0.5, 0.4)
    ax2.grid(axis="x", linestyle=":", alpha=0.3)

    fig.suptitle(
        f"LIRI-JP HCC (n={len(df)}): c5 score → HLA-C ↓ — etiology- and immune-context-independent",
        fontsize=11, y=1.03,
    )
    fig.tight_layout()
    out = OUT / "figure_liri_c5_tertile"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
