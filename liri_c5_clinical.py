"""LIRI-JP c5 score vs clinical variables (n=122 with both genotype + EGA TPM).

Tests:
  c5_score ~ age           (Pearson + Spearman)
  c5_score ~ gender        (Mann-Whitney U)
  c5_score ~ viral_status  (Kruskal-Wallis across HBV / HCV / NBNC)
  c5_score ~ T_stage       (Spearman, ordinal 1-4)

Same set of tests for log2(HLA-C TPM + 1) — the only c5 axis signal that
replicated in LIRI HCC tumor.

Stratified: HLA-C ~ c5 within each viral_status group (test whether the
c5 → HLA-C ↓ effect is etiology-dependent).

Output:
  results/liri_c5_clinical.csv
  results/figure_liri_c5_clinical.{pdf,png}
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import (
    kruskal, mannwhitneyu, pearsonr, spearmanr,
)

REPO = Path(__file__).resolve().parent
LIRI_C5 = REPO / "results/liri_c5_score.csv"     # already has c5_score + 7 genes (n=122)
CLIN = "/mnt/e/LICA-CN-analysis/LIRI-JP_analysis/data/clinical_data/EGA_clinical_matched.csv"
OUT = REPO / "results"


def main():
    sys.stdout.reconfigure(line_buffering=True)
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(LIRI_C5)
    print(f"liri c5 table: {len(df)} samples, columns: {list(df.columns)}")

    clin = pd.read_csv(CLIN)
    print(f"clinical table: {len(clin)} rows")
    clin = clin.rename(columns={"clinical_sample_id": "rk_id"})
    keep = ["rk_id", "age", "gender", "viral_status",
             "T_stage", "N_stage", "M_stage"]
    clin = clin[keep]

    merged = df.merge(clin, on="rk_id", how="inner")
    print(f"merged: {len(merged)} samples with clinical")
    print(f"  viral_status counts:\n{merged['viral_status'].value_counts()}")
    print(f"  gender counts:\n{merged['gender'].value_counts()}")
    print(f"  T_stage counts:\n{merged['T_stage'].value_counts()}")

    # Outcome variables
    merged["log2_HLA_C"] = np.log2(merged["HLA-C"].astype(float) + 1)
    merged["log2_HLA_B"] = np.log2(merged["HLA-B"].astype(float) + 1)
    merged["log2_MICA"]  = np.log2(merged["MICA"].astype(float) + 1)
    merged["log2_MICB"]  = np.log2(merged["MICB"].astype(float) + 1)

    rows = []
    outcomes = [("c5_score", merged["c5_score"].astype(float)),
                ("log2_HLA-C", merged["log2_HLA_C"]),
                ("log2_HLA-B", merged["log2_HLA_B"]),
                ("log2_MICA",  merged["log2_MICA"]),
                ("log2_MICB",  merged["log2_MICB"])]

    for ovar, y in outcomes:
        # 1. age (continuous)
        x = merged["age"].astype(float).values
        v = ~(np.isnan(x) | np.isnan(y.values))
        if v.sum() >= 5:
            rp, pp = pearsonr(x[v], y.values[v])
            rs, ps = spearmanr(x[v], y.values[v])
        else:
            rp = pp = rs = ps = float("nan")
        rows.append({"outcome": ovar, "predictor": "age",
                      "test": "Pearson/Spearman",
                      "stat": rp, "p": pp, "stat2": rs, "p2": ps,
                      "n": int(v.sum())})

        # 2. gender (M vs F)
        m = merged["gender"] == "M"
        f = merged["gender"] == "F"
        ym = y.values[m & ~y.isna().values]; yf = y.values[f & ~y.isna().values]
        if len(ym) >= 5 and len(yf) >= 5:
            u, pu = mannwhitneyu(ym, yf, alternative="two-sided")
            mean_M = float(ym.mean()); mean_F = float(yf.mean())
        else:
            u = pu = mean_M = mean_F = float("nan")
        rows.append({"outcome": ovar, "predictor": "gender(M-F)",
                      "test": "MannWhitney",
                      "stat": mean_M - mean_F, "p": pu,
                      "stat2": float("nan"), "p2": float("nan"),
                      "n": int((m | f).sum())})

        # 3. viral_status — Kruskal across HBV/HCV/NBNC
        groups = {}
        for vs in ["HBV", "HCV", "NBNC"]:
            sel = merged["viral_status"] == vs
            yv = y.values[sel.values & ~y.isna().values]
            if len(yv) >= 3:
                groups[vs] = yv
        if len(groups) >= 2:
            stat, pkw = kruskal(*groups.values())
            means = {k: float(v.mean()) for k, v in groups.items()}
            mean_str = " | ".join(f"{k}:{v:.2f}" for k, v in means.items())
        else:
            stat = pkw = float("nan"); mean_str = ""
        rows.append({"outcome": ovar, "predictor": "viral_status(HBV/HCV/NBNC)",
                      "test": "KruskalWallis",
                      "stat": stat, "p": pkw,
                      "stat2": float("nan"), "p2": float("nan"),
                      "n": sum(len(v) for v in groups.values()),
                      "means": mean_str})

        # 4. T_stage (ordinal, Spearman)
        x = merged["T_stage"].astype(float).values
        v = ~(np.isnan(x) | np.isnan(y.values))
        if v.sum() >= 5:
            rs, ps = spearmanr(x[v], y.values[v])
        else:
            rs = ps = float("nan")
        rows.append({"outcome": ovar, "predictor": "T_stage",
                      "test": "Spearman",
                      "stat": rs, "p": ps,
                      "stat2": float("nan"), "p2": float("nan"),
                      "n": int(v.sum())})

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "liri_c5_clinical.csv", index=False)
    print("\n=== LIRI clinical correlations ===")
    print(res.to_string(index=False, formatters={
        "stat": "{:+.3f}".format, "p": "{:.4f}".format,
        "stat2": "{:+.3f}".format, "p2": "{:.4f}".format,
    }))

    # Stratified: HLA-C ~ c5 within each viral_status
    print("\n=== HLA-C ~ c5 score, stratified by viral_status ===")
    for vs in ["HBV", "HCV", "NBNC"]:
        sub = merged[merged["viral_status"] == vs]
        if len(sub) < 5: continue
        x = sub["c5_score"].astype(float).values
        y = sub["log2_HLA_C"].values
        v = ~(np.isnan(x) | np.isnan(y))
        if v.sum() < 5: continue
        rp, pp = pearsonr(x[v], y[v])
        rs, ps = spearmanr(x[v], y[v])
        print(f"  {vs:5s} n={v.sum():3d}  "
              f"Pearson r={rp:+.3f} (p={pp:.3f})  "
              f"Spearman ρ={rs:+.3f} (p={ps:.3f})")

    # ===== Figure: 5 outcomes × 4 predictors panel =====
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })

    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5))
    pal_v = {"HBV":"#d65454","HCV":"#3a8c5a","NBNC":"#5a4fcf"}

    # Row 1: c5 score vs (age, gender, viral, T_stage)
    # Row 2: log2(HLA-C) vs (age, gender, viral, T_stage)
    for r_idx, (oname, yvals) in enumerate(
            [("c5 score", merged["c5_score"].astype(float)),
             ("log2(HLA-C TPM+1)", merged["log2_HLA_C"])]):
        # age
        ax = axes[r_idx, 0]
        x = merged["age"].astype(float).values
        v = ~(np.isnan(x) | yvals.isna().values)
        ax.scatter(x[v], yvals.values[v], s=22, color="#444", alpha=0.6,
                    edgecolor="#222", linewidth=0.4)
        if v.sum() >= 5:
            m, b = np.polyfit(x[v], yvals.values[v], 1)
            xs = np.linspace(x[v].min(), x[v].max(), 50)
            ax.plot(xs, m*xs + b, color="#a33", lw=0.9, ls="--")
            r, p = pearsonr(x[v], yvals.values[v])
            ax.set_title(f"{oname} ~ age  r={r:+.2f}, p={p:.3f}",
                          loc="left", fontsize=9.5)
        ax.set_xlabel("age (years)"); ax.set_ylabel(oname)

        # gender
        ax = axes[r_idx, 1]
        for i, g in enumerate(["F", "M"]):
            sel = merged["gender"] == g
            yv = yvals.values[sel.values & ~yvals.isna().values]
            ax.scatter(np.full(len(yv), i) + np.random.uniform(-0.1, 0.1, len(yv)),
                        yv, s=22, color=pal_v.get(g, "#666"), alpha=0.55,
                        edgecolor="#222", linewidth=0.4)
            ax.scatter([i], [np.median(yv)], s=80, marker="_",
                        color="#222", linewidth=2)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["F", "M"])
        ax.set_xlabel("gender")
        u_p = res[(res["outcome"] == ("c5_score" if r_idx == 0 else "log2_HLA-C"))
                  & (res["predictor"] == "gender(M-F)")]["p"].iloc[0]
        ax.set_title(f"{oname} ~ gender  MannWhitney p={u_p:.3f}",
                      loc="left", fontsize=9.5)
        ax.set_ylabel(oname)

        # viral_status
        ax = axes[r_idx, 2]
        positions = []; labels = []
        for i, vs in enumerate(["HBV", "HCV", "NBNC"]):
            sel = merged["viral_status"] == vs
            yv = yvals.values[sel.values & ~yvals.isna().values]
            if len(yv) == 0: continue
            xj = np.full(len(yv), i) + np.random.uniform(-0.12, 0.12, len(yv))
            ax.scatter(xj, yv, s=22, color=pal_v[vs], alpha=0.6,
                        edgecolor="#222", linewidth=0.4, label=f"{vs} n={len(yv)}")
            ax.scatter([i], [np.median(yv)], s=80, marker="_",
                        color="#222", linewidth=2)
            positions.append(i); labels.append(f"{vs}\nn={len(yv)}")
        ax.set_xticks(positions); ax.set_xticklabels(labels, fontsize=8)
        kw_p = res[(res["outcome"] == ("c5_score" if r_idx == 0 else "log2_HLA-C"))
                   & (res["predictor"] == "viral_status(HBV/HCV/NBNC)")]["p"].iloc[0]
        ax.set_title(f"{oname} ~ viral_status  KW p={kw_p:.3f}",
                      loc="left", fontsize=9.5)
        ax.set_xlabel("viral_status"); ax.set_ylabel(oname)

        # T_stage
        ax = axes[r_idx, 3]
        x = merged["T_stage"].astype(float).values
        v = ~(np.isnan(x) | yvals.isna().values)
        for ts in [1, 2, 3, 4]:
            sel = (x == ts) & v
            yv = yvals.values[sel]
            if len(yv) == 0: continue
            xj = np.full(len(yv), ts) + np.random.uniform(-0.12, 0.12, len(yv))
            ax.scatter(xj, yv, s=22, color="#7e3da7", alpha=0.55,
                        edgecolor="#222", linewidth=0.4)
            ax.scatter([ts], [np.median(yv)], s=80, marker="_",
                        color="#222", linewidth=2)
        rs, ps = spearmanr(x[v], yvals.values[v])
        ax.set_xticks([1, 2, 3, 4])
        ax.set_xlabel("T stage"); ax.set_ylabel(oname)
        ax.set_title(f"{oname} ~ T_stage  Spearman ρ={rs:+.2f}, p={ps:.3f}",
                      loc="left", fontsize=9.5)

    fig.suptitle(
        f"LIRI-JP c5 axis: clinical correlations (n={len(merged)} HCC samples)\n"
        "Top: c5 score (Σ alt dosage of c5 top-5% SNVs).  "
        "Bottom: log2(HLA-C TPM+1) — the c5 axis component that replicated in LIRI.",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    out = OUT / "figure_liri_c5_clinical"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
