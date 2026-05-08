"""LIRI-JP survival analysis: c5 score → outcome (Cox + KM).

Prognosis legend (Fujimoto 2016 S-table1, footnote f):
  0 = alive at last follow-up (censored)
  1 = death from liver cancer  (DSS event)
  2 = surgery-related death    (competing risk)
  3 = death from other diseases (competing risk)

Two endpoints tested:
  DSS (disease-specific survival): event = (Prognosis == 1)
  OS  (overall survival):           event = (Prognosis >= 1)

Models:
  Univariate:  c5_score (z-scored) only
  Adjusted:    + T_stage + viral_status + age + sex
  HLA-C model: + log2(HLA-C TPM+1)  (does c5 act through HLA-C?)

KM stratification: c5 score tertile (low / mid / high).

Output:
  results/liri_c5_survival_cox.csv
  results/figure_liri_c5_survival.{pdf,png}
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

REPO = Path(__file__).resolve().parent
LIRI_C5 = REPO / "results/liri_c5_score.csv"
CLIN = "/mnt/e/LICA-CN-analysis/LIRI-JP_analysis/data/clinical_data/EGA_clinical_matched.csv"
FUJIMOTO = REPO / "results/_fujimoto_stable1.csv"
OUT = REPO / "results"


def main():
    sys.stdout.reconfigure(line_buffering=True)
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(LIRI_C5)
    clin = pd.read_csv(CLIN).rename(columns={"clinical_sample_id": "rk_id"})
    df = df.merge(clin[["rk_id", "age", "gender", "viral_status",
                         "T_stage"]], on="rk_id")

    fuji = pd.read_csv(FUJIMOTO)
    fuji = fuji.rename(columns={"ID": "rk_id",
                                 "Prognosisf": "prognosis",
                                 "Overall survival (month)": "OS_month"})
    df = df.merge(fuji[["rk_id", "prognosis", "OS_month"]], on="rk_id", how="left")
    df = df[df["OS_month"].notna() & df["prognosis"].notna()].copy()
    df["OS_month"] = df["OS_month"].astype(float)
    df["prognosis"] = df["prognosis"].astype(int)
    print(f"effective n with survival: {len(df)}")
    print(f"  prognosis distribution:\n{df['prognosis'].value_counts().sort_index()}")
    print(f"  OS month: median {df['OS_month'].median():.1f},  "
           f"range [{df['OS_month'].min():.1f}, {df['OS_month'].max():.1f}]")

    # Encode events
    df["DSS_event"] = (df["prognosis"] == 1).astype(int)   # liver-cancer death only
    df["OS_event"]  = (df["prognosis"] >= 1).astype(int)   # any death
    print(f"  DSS events: {df['DSS_event'].sum()}/{len(df)}")
    print(f"  OS events:  {df['OS_event'].sum()}/{len(df)}")

    # Encode covariates
    cs = df["c5_score"].astype(float)
    df["c5_score_z"] = (cs - cs.mean()) / cs.std(ddof=0)
    df["log2_HLA_C"] = np.log2(df["HLA-C"].astype(float) + 1)
    df["sex_M"] = (df["gender"] == "M").astype(int)
    df["viral_HCV"]  = (df["viral_status"] == "HCV").astype(int)
    df["viral_NBNC"] = (df["viral_status"] == "NBNC").astype(int)
    df["age"] = df["age"].astype(float)
    df["T_stage"] = df["T_stage"].astype(float)
    df = df.dropna(subset=["c5_score_z", "T_stage", "age", "viral_status"])
    print(f"  n after covariate dropna: {len(df)}")

    # Tertile of c5_score for KM
    q33, q66 = np.quantile(df["c5_score"], [1/3, 2/3])
    df["c5_tertile"] = pd.cut(df["c5_score"], bins=[-np.inf, q33, q66, np.inf],
                                labels=["low", "mid", "high"])

    # ===== Cox regression =====
    rows = []
    cox_models = {
        "M1: c5 only":           ["c5_score_z"],
        "M2: + T_stage":         ["c5_score_z", "T_stage"],
        "M3: + T + viral":       ["c5_score_z", "T_stage", "viral_HCV", "viral_NBNC"],
        "M4: + T + viral + age + sex": ["c5_score_z", "T_stage", "viral_HCV", "viral_NBNC", "age", "sex_M"],
        "M5: M4 + log2(HLA-C)":  ["c5_score_z", "T_stage", "viral_HCV", "viral_NBNC", "age", "sex_M", "log2_HLA_C"],
    }

    for endpoint, ev_col in [("DSS", "DSS_event"), ("OS", "OS_event")]:
        for label, covs in cox_models.items():
            sub = df[["OS_month", ev_col] + covs].dropna()
            cph = CoxPHFitter(penalizer=0.001)   # tiny penalty to handle multicollinearity
            try:
                cph.fit(sub, duration_col="OS_month", event_col=ev_col)
                hr = float(np.exp(cph.params_["c5_score_z"]))
                lo = float(np.exp(cph.confidence_intervals_.loc["c5_score_z", "95% lower-bound"]))
                hi = float(np.exp(cph.confidence_intervals_.loc["c5_score_z", "95% upper-bound"]))
                p = float(cph.summary.loc["c5_score_z", "p"])
                # log2 HLA-C HR if in model
                if "log2_HLA_C" in covs:
                    hr_hlac = float(np.exp(cph.params_["log2_HLA_C"]))
                    p_hlac = float(cph.summary.loc["log2_HLA_C", "p"])
                else:
                    hr_hlac = float("nan"); p_hlac = float("nan")
            except Exception as e:
                print(f"  WARN {endpoint} {label}: {e}")
                hr = lo = hi = p = float("nan"); hr_hlac = p_hlac = float("nan")
            rows.append({
                "endpoint": endpoint, "model": label,
                "n": len(sub), "events": int(sub[ev_col].sum()),
                "c5_HR": hr, "c5_HR_lo": lo, "c5_HR_hi": hi, "c5_p": p,
                "HLA-C_HR": hr_hlac, "HLA-C_p": p_hlac,
            })

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "liri_c5_survival_cox.csv", index=False)

    print("\n=== Cox PH: c5_score (z) hazard ratio per +1 SD ===")
    print(res.to_string(index=False, formatters={
        "c5_HR": "{:.3f}".format, "c5_HR_lo": "{:.2f}".format,
        "c5_HR_hi": "{:.2f}".format, "c5_p": "{:.4f}".format,
        "HLA-C_HR": "{:.3f}".format, "HLA-C_p": "{:.4f}".format,
    }))

    # ===== KM by c5 tertile (DSS) =====
    print("\n=== KM by c5_score tertile (log-rank across 3 strata) ===")
    for endpoint, ev_col in [("DSS", "DSS_event"), ("OS", "OS_event")]:
        groups = [df[df["c5_tertile"] == t] for t in ["low", "mid", "high"]]
        result = multivariate_logrank_test(
            df["OS_month"], df["c5_tertile"], df[ev_col])
        print(f"  {endpoint} 3-group log-rank: stat={result.test_statistic:.2f}, "
              f"p={result.p_value:.4f}")
        # Pairwise low vs high
        lo = df[df["c5_tertile"] == "low"]; hi = df[df["c5_tertile"] == "high"]
        rl = logrank_test(lo["OS_month"], hi["OS_month"],
                           event_observed_A=lo[ev_col], event_observed_B=hi[ev_col])
        print(f"  {endpoint} low vs high log-rank: stat={rl.test_statistic:.2f}, "
              f"p={rl.p_value:.4f}")

    # ===== Figure =====
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    pal_t = {"low":"#5a4fcf", "mid":"#888", "high":"#d65454"}

    for ax, (endpoint, ev_col) in zip(axes, [("DSS", "DSS_event"),
                                              ("OS", "OS_event")]):
        for t in ["low", "mid", "high"]:
            sub = df[df["c5_tertile"] == t]
            kmf = KaplanMeierFitter()
            kmf.fit(sub["OS_month"], event_observed=sub[ev_col],
                     label=f"c5 {t} (n={len(sub)}, ev={int(sub[ev_col].sum())})")
            kmf.plot_survival_function(ax=ax, color=pal_t[t], lw=1.6, ci_show=False)
        # Log-rank (3-group)
        result = multivariate_logrank_test(
            df["OS_month"], df["c5_tertile"], df[ev_col])
        ax.text(0.05, 0.05, f"3-group log-rank\np = {result.p_value:.3f}",
                 transform=ax.transAxes, fontsize=9,
                 bbox=dict(boxstyle="round,pad=0.3", fc="#f5f5f5",
                            ec="#888", lw=0.5))
        ax.set_xlabel("Months from surgery")
        ax.set_ylabel(f"{endpoint} probability")
        ax.set_title(f"{endpoint} by c5_score tertile (n={len(df)} HCC)",
                      loc="left")
        ax.set_ylim(0, 1.02)
        ax.grid(linestyle=":", alpha=0.3)
        ax.legend(fontsize=7.5, loc="lower left")
    fig.suptitle(
        "LIRI-JP survival by c5_score tertile (Fujimoto 2016 S-table1 OS)\n"
        "DSS = death from liver cancer.  OS = any-cause death.",
        fontsize=10.5, y=1.04,
    )
    fig.tight_layout()
    out = OUT / "figure_liri_c5_survival"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
