"""LIRI-JP c5 proxy score validation.

Build c5 score per LIRI-JP sample from c5 top-5% SNVs (NMF k=8 26-pop)
and correlate with FPKM-UQ expression of HLA-B / HLA-C / MICA / MICB /
HCG27 / HCP5.

Predictions (from GTEx Whole_Blood + Liver):
  HLA-B  ↑   (positive correlation with c5 score)
  HLA-C  ↓   (negative)
  HCG27  ↑   (positive)
  MICA   ≈0  (no MICA eQTL signature in c5)
  MICB   ↑   (positive)

Output:
  results/liri_c5_score.csv     (per-sample score + expression)
  results/liri_c5_correlation.csv
  results/figure_liri_c5_validation.{pdf,png}
"""

from __future__ import annotations

import os

import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

REPO = Path(__file__).resolve().parent
H_PARQUET = REPO / "results/nmf_H_26_k8.parquet"
LIRI_VCF = str(Path(os.environ.get("LIRI_DATA_DIR", "data/liri_jp")) / "pcawg_germline/LIRI-JP_MHC_germline.vcf.gz")
# EGA-derived RNA-seq, gene-symbol indexed TPM, 130 RK samples (intersect 122 with genotype set).
LIRI_EXPR = str(Path(os.environ.get("LIRI_DATA_DIR", "data/liri_jp")) / "results/expression_matrix/tpm_gene_named.csv")
LIRI_EXPR_UNIT = "TPM"
# Maps VCF sample UUID -> RK ID (225 matched samples).
LIRI_VCF_TO_RK = str(Path(os.environ.get("LIRI_DATA_DIR", "data/liri_jp")) / "pcawg_germline/LIRI-JP_rs11509487_genotypes.tsv")
OUT = REPO / "results"

C5_COMPONENT = 5
NMF_TOP_FRAC = 0.05

# Gene symbols (EGA-derived expression matrix is symbol-indexed)
TARGETS = ["HLA-B", "HLA-C", "MICA", "MICB", "HCG27", "HCP5", "HLA-A"]
# HLA-A = negative control (unrelated to c5)


def main():
    sys.stdout.reconfigure(line_buffering=True)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)

    # 1. Load c5 top-5% positions (NMF k=8 26-pop)
    H = pd.read_parquet(H_PARQUET)
    snv_pos = np.array([int(c) for c in H.columns])
    h_mat = H.values
    top_n = max(1, int(NMF_TOP_FRAC * len(snv_pos)))
    c5_idx = np.argsort(h_mat[C5_COMPONENT])[-top_n:]
    c5_pos = sorted(int(p) for p in snv_pos[c5_idx])
    c5_weights = dict(zip([int(p) for p in snv_pos[c5_idx]],
                            h_mat[C5_COMPONENT][c5_idx]))
    print(f"c5 top-{top_n} SNV positions (GRCh37); range "
           f"{c5_pos[0]:,}-{c5_pos[-1]:,}")

    # 2. Subset LIRI VCF to those positions
    print(f"\nsubsetting LIRI VCF...")
    pos_str = ",".join(f"6:{p}-{p}" for p in c5_pos)
    # Use a temp BED to avoid arg overflow
    bed_path = OUT / "_c5_positions.bed"
    with open(bed_path, "w") as fh:
        for p in c5_pos:
            fh.write(f"6\t{p-1}\t{p}\n")
    import subprocess
    cmd = ["bcftools", "view", "-R", str(bed_path), "-Ou", LIRI_VCF]
    extract = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    cmd2 = ["bcftools", "query", "-f",
            "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
            "-l", "-"]
    # Two-step: first get samples, then genotypes
    samples = subprocess.check_output(["bcftools", "query", "-l", LIRI_VCF],
                                        text=True).strip().split("\n")
    print(f"  LIRI VCF: {len(samples)} samples")

    # Use bcftools query directly with -f
    proc = subprocess.run(
        ["bcftools", "query", "-R", str(bed_path),
         "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", LIRI_VCF],
        capture_output=True, text=True, timeout=300,
    )
    if proc.returncode != 0:
        print("STDERR:", proc.stderr[:500])
        raise SystemExit(proc.returncode)
    lines = proc.stdout.strip().split("\n")
    print(f"  records returned: {len(lines)}")

    # Parse genotype matrix
    sample_idx = {s: i for i, s in enumerate(samples)}
    n_samples = len(samples)
    geno_rows = []
    pos_kept = []
    n_skipped = 0
    for line in lines:
        f = line.split("\t")
        chrom = f[0]; pos = int(f[1]); ref = f[2]; alt = f[3]
        if "," in alt or alt.startswith("<"): continue
        if len(ref) != 1 or len(alt) != 1:
            n_skipped += 1; continue
        gts = f[4:]
        if len(gts) != n_samples:
            n_skipped += 1; continue
        # Convert to alt-allele dosage (0/1/2 or NaN)
        dose = np.full(n_samples, np.nan)
        for i, g in enumerate(gts):
            if "|" in g: a, b = g.split("|", 1)
            elif "/" in g: a, b = g.split("/", 1)
            else: continue
            if a == "." or b == ".": continue
            try: dose[i] = int(a) + int(b)
            except ValueError: continue
        geno_rows.append((pos, dose))
        pos_kept.append(pos)
    print(f"  c5 SNVs found in LIRI: {len(pos_kept)} / {len(c5_pos)} requested")
    print(f"  (skipped: {n_skipped} — multiallelic / indel)")

    if len(pos_kept) == 0:
        raise SystemExit("no c5 SNVs in LIRI VCF")

    # Build sample × SNV dosage matrix
    geno_mat = np.array([row[1] for row in geno_rows]).T   # (n_samples, n_snv)
    print(f"  genotype matrix shape: {geno_mat.shape}")
    print(f"  per-SNV missingness: median={np.median(np.isnan(geno_mat).mean(axis=0)):.3f}")

    # 3. c5 score = sum of alt dosage (binary count style)
    # Skip SNVs with high missingness (>20%); drop samples with high missing too
    snv_miss = np.isnan(geno_mat).mean(axis=0)
    keep_snv = snv_miss <= 0.2
    geno = geno_mat[:, keep_snv]
    pos_use = np.array(pos_kept)[keep_snv]
    print(f"  SNVs after missingness filter: {keep_snv.sum()}")

    # Imputation with SNV-mean for the remaining missing genotypes
    snv_mean = np.nanmean(geno, axis=0)
    inds_nan = np.where(np.isnan(geno))
    geno[inds_nan] = snv_mean[inds_nan[1]]

    c5_score = geno.sum(axis=1)
    sample_score = pd.DataFrame({"uuid": samples, "c5_score": c5_score})

    # 4. VCF UUID → RK mapping (from rs11509487 genotype TSV)
    umap = pd.read_csv(LIRI_VCF_TO_RK, sep="\t",
                        usecols=["liriju_sample_id", "vcf_sample_id"])
    umap = umap.rename(columns={"vcf_sample_id": "uuid",
                                 "liriju_sample_id": "rk_id"})
    sample_score = sample_score.merge(umap, on="uuid", how="left")
    print(f"\n  samples with RK ID: {sample_score['rk_id'].notna().sum()}")

    # 5. Load expression — EGA-derived matrix uses gene symbols as row index;
    # first column is unlabeled (gene name), rest are RK sample IDs.
    print(f"\nloading expression matrix ({LIRI_EXPR_UNIT})...")
    expr = pd.read_csv(LIRI_EXPR)
    gene_col = expr.columns[0]   # often unnamed -> 'Unnamed: 0' or empty
    expr_samples = [c for c in expr.columns[1:]]
    print(f"  expression samples: {len(expr_samples)}")
    feats = expr[gene_col].astype(str)

    target_rows = {}
    for gene_sym in TARGETS:
        m = (feats == gene_sym)
        if m.any():
            target_rows[gene_sym] = expr.loc[m].iloc[0]
        else:
            print(f"  WARN: {gene_sym} not found in expression matrix")

    # 6. Intersect samples
    sample_score_with_rk = sample_score.dropna(subset=["rk_id"])
    expr_intersect = [s for s in expr_samples if s in sample_score_with_rk["rk_id"].values]
    print(f"\n  intersected samples (geno + expr): {len(expr_intersect)}")

    if len(expr_intersect) < 5:
        print("  TOO FEW INTERSECT — abort")
        return

    # Build merged dataframe
    score_lookup = dict(zip(sample_score_with_rk["rk_id"],
                              sample_score_with_rk["c5_score"]))
    df = pd.DataFrame({"rk_id": expr_intersect})
    df["c5_score"] = df["rk_id"].map(score_lookup)
    for gene_sym in TARGETS:
        if gene_sym in target_rows:
            row = target_rows[gene_sym]
            df[gene_sym] = [row[s] for s in expr_intersect]
        else:
            df[gene_sym] = np.nan
    df.to_csv(OUT / "liri_c5_score.csv", index=False)
    print(f"  wrote {OUT / 'liri_c5_score.csv'}")

    # 7. Correlate
    print(f"\n=== c5 score vs gene expression (n={len(df)}) ===")
    rows = []
    for gene_sym in TARGETS:
        if gene_sym not in df.columns: continue
        x = df["c5_score"].astype(float).values
        y = df[gene_sym].astype(float).values
        valid = ~(np.isnan(x) | np.isnan(y))
        if valid.sum() < 5: continue
        x_v = x[valid]; y_v = y[valid]
        # log-transform expression (TPM; +1 to handle zeros)
        y_log = np.log2(y_v + 1)
        try:
            r_p, p_p = pearsonr(x_v, y_log)
            r_s, p_s = spearmanr(x_v, y_log)
        except Exception:
            r_p = p_p = r_s = p_s = float("nan")
        # Median expression high vs low c5 (median split)
        med = np.median(x_v)
        high = y_log[x_v >= med]; low = y_log[x_v < med]
        rows.append({
            "gene": gene_sym,
            "n_samples": int(valid.sum()),
            "pearson_r_log2_FPKM": r_p, "pearson_p": p_p,
            "spearman_rho": r_s, "spearman_p": p_s,
            f"log2_{LIRI_EXPR_UNIT}_high": float(high.mean()) if len(high) else float("nan"),
            f"log2_{LIRI_EXPR_UNIT}_low": float(low.mean()) if len(low) else float("nan"),
            "diff_high_minus_low": float(high.mean() - low.mean()) if len(high) and len(low) else float("nan"),
        })
    cor = pd.DataFrame(rows)
    cor.to_csv(OUT / "liri_c5_correlation.csv", index=False)
    print(cor.to_string(index=False, formatters={
        "pearson_r_log2_FPKM": "{:+.3f}".format, "pearson_p": "{:.3f}".format,
        "spearman_rho": "{:+.3f}".format, "spearman_p": "{:.3f}".format,
        f"log2_{LIRI_EXPR_UNIT}_high": "{:.2f}".format,
        f"log2_{LIRI_EXPR_UNIT}_low": "{:.2f}".format,
        "diff_high_minus_low": "{:+.3f}".format,
    }))

    # ===== Figure: scatter c5 score vs each gene =====
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    palette = {"HLA-B":"#5a4fcf","HLA-C":"#3a8c5a","MICA":"#d65454",
                "MICB":"#b5722d","HCG27":"#7e3da7","HCP5":"#888",
                "HLA-A":"#1f4e8c"}
    n_genes = sum(1 for g in TARGETS if g in df.columns)
    fig, axes = plt.subplots(2, 4, figsize=(15, 7))
    for ax, gene_sym in zip(axes.flat, list(TARGETS) + [None]*max(0, 8 - len(TARGETS))):
        if gene_sym is None or gene_sym not in df.columns:
            ax.axis("off"); continue
        x = df["c5_score"].astype(float).values
        y = df[gene_sym].astype(float).values
        valid = ~(np.isnan(x) | np.isnan(y))
        if valid.sum() < 5:
            ax.text(0.5, 0.5, f"{gene_sym}\n(n<5)", transform=ax.transAxes,
                    ha="center", va="center"); ax.axis("off"); continue
        y_log = np.log2(y[valid] + 1)
        ax.scatter(x[valid], y_log, s=22, color=palette.get(gene_sym, "#444"),
                    alpha=0.65, edgecolor="#222", linewidth=0.4)
        # fit line
        m, b = np.polyfit(x[valid], y_log, 1)
        xs = np.linspace(x[valid].min(), x[valid].max(), 50)
        ax.plot(xs, m*xs + b, color="#222", lw=0.8, ls="--")
        # stats
        try:
            r_p, p_p = pearsonr(x[valid], y_log)
        except Exception:
            r_p = p_p = float("nan")
        ax.set_title(f"{gene_sym}  r={r_p:+.2f}, p={p_p:.3g}",
                      loc="left", fontsize=10)
        ax.set_xlabel("c5 score (Σ alt dosage)")
        ax.set_ylabel(rf"$\log_2$({LIRI_EXPR_UNIT} + 1)")
        ax.grid(linestyle=":", alpha=0.3)
    fig.suptitle(
        f"LIRI-JP: c5 proxy score vs HLA region expression (n={len(df)} HCC samples)\n"
        f"Predictions: HLA-B↑, HLA-C↓, HCG27↑, MICA≈0, MICB↑",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    out = OUT / "figure_liri_c5_validation"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)
    print(f"\nTotal: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
