#!/usr/bin/env python3
"""Robustness Task 1 — LD-clumping to independent-locus signal counts.

Pseudoreplication check: the NMF component eQTL enrichment / sign-coherence is
computed over top-5% signature SNVs that may be in LD. Here we greedily LD-clump
each component's top-5% SNVs (by descending H-loading) and re-evaluate eQTL
enrichment and sign-coherence at (A) index-only and (B) clump level, to confirm
the c5 = HLA-B↑/HLA-C↓ and c4/c6 = MICA↓ signals survive at independent loci.

This is a STANDALONE robustness/sensitivity analysis. It does NOT modify or
regenerate any existing figure, table, notebook, or result file; it only writes
new results/robustness_ldclump_*.tsv files.

Inputs (GRCh37; defaults are repo-relative, override with argparse):
  results/nmf_H_26_k8.parquet            (8 x 7116 H-loadings; cols = b37 pos)
  results/_supp_hap_matrix_26pop.npy     (2111 carrier haplotypes x 7116 SNV)
  results/_supp_hap_positions.npy        (7116 b37 positions, aligned to H cols)
  results/_supp_hap_meta.parquet         (per-haplotype pop / super_pop)
  results/eqtl_targets_per_gene.csv      (gene, b37_pos, tissue) GTEx eQTL hits
  results/_axis_nes_full.csv             (tissue, gene, component, b37_pos, NES)

Outputs:
  results/robustness_ldclump_component_summary.tsv
  results/robustness_ldclump_eqtl_enrichment.tsv
  results/robustness_ldclump_sign_coherence.tsv
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SEED = 0
TOP_FRAC = 0.05
R2_PRIMARY = 0.2
R2_SENS = 0.1
NES_TISSUE = "Whole_Blood"          # primary tissue for direction (matches panel a/b)

# component -> target genes to evaluate; c1 is a comparison control
COMP_GENES = {
    "c5": ["HLA-B", "HLA-C"],
    "c4": ["MICA"],
    "c6": ["MICA"],
    "c1": ["HLA-B", "HLA-C", "MICA"],   # control
}
EXPECTED_DIR = {("c5", "HLA-B"): +1, ("c5", "HLA-C"): -1,
                ("c4", "MICA"): -1, ("c6", "MICA"): -1}


def greedy_clump(order_idx, r2mat, thr):
    """order_idx: SNV indices sorted by descending H-loading (into r2mat space).
    r2mat: square r^2 matrix over those SNVs. Returns list of clumps; each clump
    is a list of original-order positions, clump[0] is the index SNV."""
    n = len(order_idx)
    used = np.zeros(n, dtype=bool)
    clumps = []
    for i in range(n):
        if used[i]:
            continue
        members = [i]
        used[i] = True
        for j in range(n):
            if not used[j] and r2mat[i, j] >= thr:
                used[j] = True
                members.append(j)
        clumps.append([order_idx[m] for m in members])
    return clumps


def r2_matrix(X):
    """X: (n_hap, n_snv) 0/1 matrix -> (n_snv, n_snv) r^2 (pearson^2).
    Constant columns -> r^2 = 0 against everything (no LD)."""
    Xc = X - X.mean(axis=0, keepdims=True)
    sd = Xc.std(axis=0)
    nz = sd > 0
    R = np.zeros((X.shape[1], X.shape[1]), dtype=float)
    if nz.sum() >= 2:
        sub = Xc[:, nz] / sd[nz]
        cc = (sub.T @ sub) / X.shape[0]
        idx = np.where(nz)[0]
        R[np.ix_(idx, idx)] = cc ** 2
    np.fill_diagonal(R, 1.0)
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--out-dir", default="results")
    args = ap.parse_args()
    np.random.seed(SEED)
    RES = Path(args.results_dir)
    OUT = Path(args.out_dir)
    OUT.mkdir(parents=True, exist_ok=True)

    # ---- load ----
    H = pd.read_parquet(RES / "nmf_H_26_k8.parquet")
    h_pos = np.array([int(c) for c in H.columns])
    h_mat = H.values
    K, n_snv = h_mat.shape
    top_n = max(1, int(TOP_FRAC * n_snv))

    hap = np.load(RES / "_supp_hap_matrix_26pop.npy").astype(np.float32)
    hap_pos = np.array([int(p) for p in np.load(RES / "_supp_hap_positions.npy",
                                                allow_pickle=True)])
    assert np.array_equal(hap_pos, h_pos), "hap positions must align to H columns"
    meta = pd.read_parquet(RES / "_supp_hap_meta.parquet")
    super_pop = meta["super_pop"].values

    eqtl = pd.read_csv(RES / "eqtl_targets_per_gene.csv")
    eqtl_sets = {g: set(eqtl[eqtl["gene"] == g]["b37_pos"].astype(int))
                 for g in eqtl["gene"].unique()}
    # window-level eQTL density per gene (background for binomial enrichment)
    bg_rate = {g: len(set(eqtl_sets[g]) & set(h_pos)) / n_snv for g in eqtl_sets}

    nes = pd.read_csv(RES / "_axis_nes_full.csv")
    nes = nes[nes["tissue"] == NES_TISSUE]
    nes_map = {}   # (component, gene, b37_pos) -> NES
    for _, r in nes.iterrows():
        nes_map[(r["component"], r["gene"], int(r["b37_pos"]))] = float(r["NES"])

    pos_index = {int(p): i for i, p in enumerate(h_pos)}

    def clump_on(rows, comp, thr):
        """LD-clump comp's top-5% SNVs using haplotype subset `rows`."""
        c = int(comp[1:])
        h = h_mat[c]
        top_order = np.argsort(h)[::-1][:top_n]            # global col idx, desc H
        sub_pos = h_pos[top_order]
        X = hap[np.ix_(rows, top_order)]
        R = r2_matrix(X)
        clumps_local = greedy_clump(list(range(len(top_order))), R, thr)
        # map back to global positions
        return [[int(sub_pos[m]) for m in cl] for cl in clumps_local], list(sub_pos)

    all_rows = np.arange(hap.shape[0])
    eas_rows = np.where(super_pop == "EAS")[0]
    eur_rows = np.where(super_pop == "EUR")[0]

    # ============ component summary (counts) ============
    summ = []
    clumps_primary = {}    # comp -> clumps (pooled, r2=0.2)
    for comp in COMP_GENES:
        cl02, top_pos = clump_on(all_rows, comp, R2_PRIMARY)
        cl01, _ = clump_on(all_rows, comp, R2_SENS)
        cl02_eas, _ = clump_on(eas_rows, comp, R2_PRIMARY)
        cl02_eur, _ = clump_on(eur_rows, comp, R2_PRIMARY)
        clumps_primary[comp] = cl02
        summ.append({
            "component": comp,
            "raw_top5pct_snv": len(top_pos),
            "n_clumps_r2_0.2": len(cl02),
            "n_clumps_r2_0.1": len(cl01),
            "n_clumps_r2_0.2_EAS": len(cl02_eas),
            "n_clumps_r2_0.2_EUR": len(cl02_eur),
        })
    summ = pd.DataFrame(summ)
    summ.to_csv(OUT / "robustness_ldclump_component_summary.tsv", sep="\t", index=False)

    # ============ enrichment + sign at index-only and clump-level ============
    enr_rows, sign_rows = [], []
    for comp, genes in COMP_GENES.items():
        clumps = clumps_primary[comp]
        n_clumps = len(clumps)
        index_snvs = [cl[0] for cl in clumps]
        for gene in genes:
            gset = eqtl_sets.get(gene, set())
            # ---- enrichment ----
            idx_overlap = sum(1 for p in index_snvs if p in gset)
            clump_overlap = sum(1 for cl in clumps if any(p in gset for p in cl))
            p_bg = bg_rate.get(gene, 0.0)
            # index-level binomial vs window background eQTL density
            binom_p = stats.binomtest(idx_overlap, n_clumps, p_bg,
                                      alternative="greater").pvalue if n_clumps else float("nan")
            enr_rows.append({
                "component": comp, "gene": gene,
                "n_clumps": n_clumps,
                "index_overlap": idx_overlap,
                "index_overlap_rate": round(idx_overlap / n_clumps, 4) if n_clumps else float("nan"),
                "clump_overlap": clump_overlap,
                "clump_overlap_rate": round(clump_overlap / n_clumps, 4) if n_clumps else float("nan"),
                "window_bg_eqtl_rate": round(p_bg, 4),
                "index_binom_p_greater": binom_p,
            })
            # ---- sign coherence ----
            exp = EXPECTED_DIR.get((comp, gene), None)
            # index-only: sign of index SNV's own NES (if it overlaps & has NES)
            idx_signs = []
            for p in index_snvs:
                v = nes_map.get((comp, gene, p))
                if v is not None and v != 0:
                    idx_signs.append(np.sign(v))
            # clump-level: per overlapping clump, sign of strongest |NES| member
            cl_signs = []
            for cl in clumps:
                cand = [(abs(nes_map[(comp, gene, p)]), np.sign(nes_map[(comp, gene, p)]))
                        for p in cl if (comp, gene, p) in nes_map and nes_map[(comp, gene, p)] != 0]
                if cand:
                    cl_signs.append(max(cand)[1])

            def coh(signs):
                if not signs:
                    return (float("nan"), 0, 0, float("nan"))
                a = np.array(signs)
                npos, nneg = int((a > 0).sum()), int((a < 0).sum())
                coherence = abs(a.mean())
                return (coherence, npos, nneg, np.sign(a.sum()))

            ic, inp, inn, idir = coh(idx_signs)
            cc, cnp, cnn, cdir = coh(cl_signs)
            sign_rows.append({
                "component": comp, "gene": gene,
                "expected_dir": exp if exp is not None else "",
                "index_n": inp + inn, "index_n_pos": inp, "index_n_neg": inn,
                "index_sign_coherence": round(ic, 3) if ic == ic else float("nan"),
                "index_dir_matches_expected": (idir == exp) if (exp and (inp + inn)) else "",
                "clump_n": cnp + cnn, "clump_n_pos": cnp, "clump_n_neg": cnn,
                "clump_sign_coherence": round(cc, 3) if cc == cc else float("nan"),
                "clump_dir_matches_expected": (cdir == exp) if (exp and (cnp + cnn)) else "",
            })

    pd.DataFrame(enr_rows).to_csv(OUT / "robustness_ldclump_eqtl_enrichment.tsv",
                                  sep="\t", index=False)
    pd.DataFrame(sign_rows).to_csv(OUT / "robustness_ldclump_sign_coherence.tsv",
                                   sep="\t", index=False)

    print("=== component summary ===")
    print(summ.to_string(index=False))
    print("\n=== enrichment (index + clump) ===")
    print(pd.DataFrame(enr_rows).to_string(index=False))
    print("\n=== sign coherence ===")
    print(pd.DataFrame(sign_rows).to_string(index=False))
    print(f"\nNES tissue = {NES_TISSUE}; r2 primary = {R2_PRIMARY}, sens = {R2_SENS}; "
          f"clumping = greedy by descending H-loading on {hap.shape[0]} carrier haplotypes.")
    print("NOTE: c1 NES not in the axis NES dump (computed for c4/c5/c6 only); "
          "c1 is reported as a count/overlap control, sign-coherence NaN by design.")


if __name__ == "__main__":
    main()
