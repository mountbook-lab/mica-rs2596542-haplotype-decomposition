"""LIRI-JP multivariate regression: c5 → HLA-C ↓ adjusted for T stage / viral / age / sex.

Reviewer concern: HLA-C correlates positively with T stage in LIRI
(Spearman ρ=+0.20, p=0.024). If higher-stage tumors over-represent
low-c5 carriers, the apparent c5 → HLA-C ↓ effect could be confounded.

Models (log2(HLA-C TPM + 1) as outcome):
  M1: c5_score  (univariate, baseline)
  M2: c5_score + T_stage
  M3: c5_score + T_stage + viral_status
  M4: c5_score + T_stage + viral_status + age + sex   (full model)

Same set for log2(HLA-B), log2(MICA), log2(MICB), log2(HCG27) for comparison.

Output:
  results/liri_c5_multivariate.csv
"""

from __future__ import annotations

import os

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent
LIRI_C5 = REPO / "results/liri_c5_score.csv"
CLIN = str(Path(os.environ.get("LIRI_DATA_DIR", "data/liri_jp")) / "data/clinical_data/EGA_clinical_matched.csv")
OUT = REPO / "results"

OUTCOMES = ["HLA-C", "HLA-B", "MICA", "MICB", "HCG27"]


def fit(Y, X, label):
    X1 = sm.add_constant(X)
    m = sm.OLS(Y, X1, missing="drop").fit()
    return {
        "model": label,
        "n": int(m.nobs),
        "c5_beta": m.params.get("c5_score_z", float("nan")),
        "c5_se":   m.bse.get("c5_score_z",   float("nan")),
        "c5_p":    m.pvalues.get("c5_score_z", float("nan")),
        "c5_t":    m.tvalues.get("c5_score_z", float("nan")),
        "r2":      float(m.rsquared),
    }


def main():
    sys.stdout.reconfigure(line_buffering=True)
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(LIRI_C5)
    clin = pd.read_csv(CLIN).rename(columns={"clinical_sample_id": "rk_id"})
    df = df.merge(clin[["rk_id", "age", "gender", "viral_status",
                         "T_stage", "N_stage", "M_stage"]], on="rk_id")
    df = df[df["T_stage"].notna() & df["age"].notna()
             & df["viral_status"].isin(["HBV", "HCV", "NBNC"])].copy()
    print(f"effective n: {len(df)}")
    print(f"  viral_status:\n{df['viral_status'].value_counts()}")

    # Standardize c5_score for interpretable coefficients
    cs = df["c5_score"].astype(float)
    df["c5_score_z"] = (cs - cs.mean()) / cs.std(ddof=0)

    # Encode categoricals
    df["sex_M"] = (df["gender"] == "M").astype(int)
    df["viral_HCV"]  = (df["viral_status"] == "HCV").astype(int)
    df["viral_NBNC"] = (df["viral_status"] == "NBNC").astype(int)
    df["age"]     = df["age"].astype(float)
    df["T_stage"] = df["T_stage"].astype(float)

    rows = []
    for gene in OUTCOMES:
        Y = np.log2(df[gene].astype(float) + 1)
        # M1 univariate
        rows.append({"outcome": gene, **fit(Y, df[["c5_score_z"]], "M1: c5 only")})
        # M2 + T_stage
        rows.append({"outcome": gene, **fit(Y,
            df[["c5_score_z", "T_stage"]], "M2: + T_stage")})
        # M3 + viral
        rows.append({"outcome": gene, **fit(Y,
            df[["c5_score_z", "T_stage", "viral_HCV", "viral_NBNC"]],
            "M3: + viral_status")})
        # M4 + age + sex (full)
        rows.append({"outcome": gene, **fit(Y,
            df[["c5_score_z", "T_stage", "viral_HCV", "viral_NBNC",
                 "age", "sex_M"]], "M4: full")})

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "liri_c5_multivariate.csv", index=False)

    print("\n=== c5_score (z-scored) coefficient on log2(gene TPM+1) ===")
    print("β = change in log2(TPM+1) per +1 SD of c5 score")
    print(res.to_string(index=False, formatters={
        "c5_beta": "{:+.3f}".format, "c5_se": "{:.3f}".format,
        "c5_t": "{:+.2f}".format, "c5_p": "{:.4f}".format,
        "r2": "{:.3f}".format,
    }))

    # HLA-C only — full model summary
    print("\n=== HLA-C full model (M4) coefficients ===")
    Y = np.log2(df["HLA-C"].astype(float) + 1)
    X = df[["c5_score_z", "T_stage", "viral_HCV", "viral_NBNC", "age", "sex_M"]]
    X1 = sm.add_constant(X)
    m = sm.OLS(Y, X1, missing="drop").fit()
    print(m.summary().tables[1])
    print(f"\nR² = {m.rsquared:.3f},  adj R² = {m.rsquared_adj:.3f},  n = {int(m.nobs)}")


if __name__ == "__main__":
    main()
