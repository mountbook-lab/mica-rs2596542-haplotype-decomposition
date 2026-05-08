"""Figure 5 — LIRI-JP HCC validation of c5 axis.

Two sub-panels labeled (a) / (b) for manuscript Figure 5:

  (a) c5 score tertile (low / mid / high, n = 39 each) vs log2(HLA-C TPM+1).
      Etiology color (HBV / HCV / NBNC). Stat box: Jonckheere trend,
      low-vs-high MWU, Spearman rho on continuous c5.

  (b) Forest plot of c5 coefficient on log2(HLA-C TPM+1) under
      progressively stricter adjustment:
        M1: + T_stage + viral_status + age + sex   (baseline)
        M2: + PTPRC                                  (general immune)
        M3: + immune_PC1 (6 markers)                 (composite immune)
        M4: + log2(HLA-A) [strict]                   (MHC class-I co-reg)

Input:
  results/liri_c5_score.csv  (n=122 with EGA TPM)
  /mnt/e/.../EGA_clinical_matched.csv
  /mnt/e/.../tpm_gene_named.csv  (immune marker TPM)

Output:
  results/figure_5_liri_panel.{pdf,png}
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


def jonckheere(values_by_group, n_perm=5000, seed=0):
    k = len(values_by_group)
    J = 0
    for i in range(k):
        for j in range(i + 1, k):
            xi = values_by_group[i]; xj = values_by_group[j]
            J += float(((xj[:, None] > xi[None, :]).sum() +
                       0.5 * (xj[:, None] == xi[None, :]).sum()))
    rng = np.random.default_rng(seed)
    all_v = np.concatenate(values_by_group)
    sizes = [len(v) for v in values_by_group]
    null = np.zeros(n_perm)
    for p in range(n_perm):
        perm = rng.permutation(all_v)
        groups = []; s = 0
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

    # Pull immune markers
    expr = pd.read_csv(EXPR)
    gene_col = expr.columns[0]
    expr_samples = list(expr.columns[1:])
    feats = expr[gene_col].astype(str)
    for g in IMMUNE_GENES:
        m = (feats == g)
        if not m.any(): continue
        row = expr.loc[m].iloc[0]
        df[g] = df["rk_id"].map(lambda rk: row[rk] if rk in expr_samples else np.nan)
        df[f"log2_{g}"] = np.log2(df[g].astype(float) + 1)
    df["log2_HLA_C"] = np.log2(df["HLA-C"].astype(float) + 1)
    df["log2_HLA_A"] = np.log2(df["HLA-A"].astype(float) + 1)

    cs = df["c5_score"].astype(float)
    df["c5_score_z"] = (cs - cs.mean()) / cs.std(ddof=0)
    q33, q66 = np.quantile(cs, [1/3, 2/3])
    df["c5_tertile"] = pd.cut(cs, bins=[-np.inf, q33, q66, np.inf],
                                labels=["low", "mid", "high"])

    # Tests
    by_t = []
    for t in ["low", "mid", "high"]:
        v = df.loc[df["c5_tertile"] == t, "log2_HLA_C"].dropna().values
        by_t.append(v)
    kw_stat, kw_p = kruskal(*by_t)
    J, p_jt = jonckheere(by_t, n_perm=5000)
    u, p_lh = mannwhitneyu(by_t[0], by_t[2], alternative="two-sided")
    rs, ps = spearmanr(cs, df["log2_HLA_C"])

    # immune PC1
    Z = df[[f"log2_{g}" for g in IMMUNE_GENES]].copy()
    for c in Z.columns:
        Z[c] = (Z[c] - Z[c].mean()) / Z[c].std(ddof=0)
    Zv = Z.dropna().values
    U, S, Vt = np.linalg.svd(Zv, full_matrices=False)
    df["immune_PC1"] = (Z.values @ Vt[0])

    df["sex_M"] = (df["gender"] == "M").astype(int)
    df["viral_HCV"]  = (df["viral_status"] == "HCV").astype(int)
    df["viral_NBNC"] = (df["viral_status"] == "NBNC").astype(int)

    Y = df["log2_HLA_C"]
    models = {
        "M1: + T + viral + age + sex":      ["c5_score_z", "T_stage", "viral_HCV",
                                              "viral_NBNC", "age", "sex_M"],
        "M2: + PTPRC":                       ["c5_score_z", "T_stage", "viral_HCV",
                                              "viral_NBNC", "age", "sex_M",
                                              "log2_PTPRC"],
        "M3: + immune_PC1":                  ["c5_score_z", "T_stage", "viral_HCV",
                                              "viral_NBNC", "age", "sex_M",
                                              "immune_PC1"],
        "M4: + log2(HLA-A) [strict]":        ["c5_score_z", "T_stage", "viral_HCV",
                                              "viral_NBNC", "age", "sex_M",
                                              "log2_HLA_A"],
    }
    rows = []
    for label, covs in models.items():
        sub = df[[Y.name] + covs].dropna()
        X1 = sm.add_constant(sub[covs])
        m = sm.OLS(sub[Y.name], X1).fit()
        rows.append({"model": label,
                      "c5_beta": float(m.params["c5_score_z"]),
                      "c5_se":   float(m.bse["c5_score_z"]),
                      "c5_p":    float(m.pvalues["c5_score_z"])})
    res = pd.DataFrame(rows)

    # ===== Figure =====
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    fig = plt.figure(figsize=(12, 5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.2])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    pal_v = {"HBV":"#d65454", "HCV":"#3a8c5a", "NBNC":"#5a4fcf"}
    rng = np.random.default_rng(0)

    # ---- (a) Tertile boxplot ----
    axA.boxplot(by_t, positions=[0, 1, 2], widths=0.5, patch_artist=False,
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
            axA.scatter(xj, ssub["log2_HLA_C"], s=22, color=pal_v[vs],
                         alpha=0.7, edgecolor="#222", linewidth=0.4,
                         label=vs if i == 0 else None)
    axA.set_xticks([0, 1, 2])
    axA.set_xticklabels([f"low\nn={len(by_t[0])}",
                          f"mid\nn={len(by_t[1])}",
                          f"high\nn={len(by_t[2])}"])
    axA.set_xlabel("c5 score tertile")
    axA.set_ylabel(r"$\log_2$(HLA-C TPM + 1)")
    axA.text(0.02, 0.98,
              f"Jonckheere-Terpstra trend p = {p_jt:.3f}\n"
              f"low vs high MWU p = {p_lh:.3f}\n"
              rf"Spearman $\rho$ (continuous c5) = {rs:+.2f}, p = {ps:.3f}",
              transform=axA.transAxes, fontsize=8.5, va="top",
              bbox=dict(boxstyle="round,pad=0.35", fc="#f8f8f8",
                         ec="#aaa", lw=0.5))
    axA.legend(title="viral status", fontsize=7.5, title_fontsize=8,
                 loc="lower left", frameon=False)
    axA.set_title("(a)  HLA-C expression by c5 score tertile",
                    loc="left", fontsize=10.5)
    axA.grid(axis="y", linestyle=":", alpha=0.3)

    # ---- (b) Forest plot ----
    labels = list(reversed([r["model"] for _, r in res.iterrows()]))
    betas  = list(reversed([r["c5_beta"] for _, r in res.iterrows()]))
    ses    = list(reversed([r["c5_se"]   for _, r in res.iterrows()]))
    ps_    = list(reversed([r["c5_p"]    for _, r in res.iterrows()]))
    ys = np.arange(len(labels))
    for y, b, s, p in zip(ys, betas, ses, ps_):
        lo, hi = b - 1.96*s, b + 1.96*s
        col = "#a33" if p < 0.05 else "#888"
        axB.errorbar([b], [y], xerr=[[b - lo], [hi - b]], fmt="o",
                       color=col, ecolor=col, capsize=4, ms=8, lw=1.5,
                       markeredgecolor="#222", markeredgewidth=0.6)
        axB.text(hi + 0.02, y,
                  rf"  $\beta$={b:+.3f}  p={p:.3f}",
                  va="center", fontsize=8.5,
                  color="#222" if p < 0.05 else "#444")
    axB.axvline(0, color="#222", lw=0.6, ls="-")
    axB.set_yticks(ys); axB.set_yticklabels(labels, fontsize=8.5)
    axB.set_xlabel(r"c5_score (z) $\beta$ on $\log_2$(HLA-C TPM+1) [95% CI]")
    axB.set_title("(b)  c5 → HLA-C robust to immune & MHC class-I co-regulation",
                    loc="left", fontsize=10.5)
    axB.set_xlim(-0.5, 0.4)
    axB.grid(axis="x", linestyle=":", alpha=0.3)

    fig.suptitle(
        f"Figure 5.  LIRI-JP HCC (n={len(df)}) validates the c5 → HLA-C ↓ axis "
        f"as etiology- and immune-context-independent",
        fontsize=11, y=1.03,
    )
    fig.tight_layout()
    out = OUT / "figure_5_liri_panel"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"wrote {out}.pdf / .png")
    plt.close(fig)


if __name__ == "__main__":
    main()
