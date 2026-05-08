"""(f) Ground-truth simulation: does NMF recover known latent branches?

Generate synthetic anchor-carrier × SNV matrix with KNOWN structure:

  4 latent branches:
    H1: AFR-enriched (60% in AFR, 10% elsewhere)
    H2: EUR-enriched (50% in EUR, 10% elsewhere)
    H3: EAS-enriched (40% in EAS, 10% elsewhere)
    H4: shared backbone (30% uniform across all populations)

  SNV pool (n=7000):
    300 SNVs are H1 markers (high frequency in H1, near-zero in other branches)
    300 SNVs are H2 markers
    300 SNVs are H3 markers
    300 SNVs are H4 (backbone) markers — present in any branch
    5800 SNVs are "noise" — random independent of branch assignment

  Hap matrix (~2000 carriers, super-pop-balanced like real data):
    each hap belongs to one branch (assigned with population-conditional prob)
    branch markers: P(carrier | branch_match) = 0.85; P(carrier | mismatch) = 0.05
    backbone: P(carrier) = 0.30 uniform
    noise: P(carrier) = MAF (random uniform 0.05-0.45)

Tests:
  1. Run UMAP+Leiden+NMF; verify recovery of H1, H2, H3, H4 components.
  2. For each NMF component, identify its top-5% SNV set and overlap with
     H1/H2/H3/H4 marker sets.
  3. Negative control: shuffle population labels (preserving frequencies);
     verify no spurious "branch" axis emerges.

Output:
  results/sim_ground_truth_summary.csv
  results/figure_simulation_ground_truth.{pdf,png}
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import igraph as ig
import leidenalg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
from scipy.optimize import linear_sum_assignment
from scipy.stats import hypergeom
from sklearn.decomposition import NMF
from sklearn.neighbors import NearestNeighbors

REPO = Path(__file__).resolve().parent
OUT = REPO / "results"

# ---- Simulation parameters ----
N_SUPER_POPS = 5  # AFR, EUR, EAS, SAS, AMR (mirroring real data)
N_PER_POP = 400   # haplotypes per super-pop  (total 2000)
N_SNV = 7000
N_MARKER_PER_BRANCH = 300
N_BACKBONE = 300

# Branch frequencies per super-pop (rows = branch, cols = super_pop)
BRANCH_FREQ = np.array([
    # AFR  EUR  EAS  SAS  AMR
    [0.60, 0.10, 0.10, 0.10, 0.10],   # H1: AFR-enriched
    [0.10, 0.50, 0.10, 0.10, 0.20],   # H2: EUR-enriched
    [0.10, 0.10, 0.40, 0.30, 0.10],   # H3: EAS/SAS-enriched
    [0.20, 0.30, 0.40, 0.50, 0.60],   # H4: backbone (residual; complementary to others)
])
# Normalise per super-pop column to sum to 1
BRANCH_FREQ = BRANCH_FREQ / BRANCH_FREQ.sum(axis=0, keepdims=True)
N_BRANCHES = BRANCH_FREQ.shape[0]

P_CARRIER_MATCH    = 0.85
P_CARRIER_MISMATCH = 0.05
P_BACKBONE_CARRIER = 0.30
NMF_TOP_FRAC = 0.05


def simulate(rng_seed: int):
    rng = np.random.default_rng(rng_seed)
    # Build super-pop labels
    super_pops = []
    for s, name in enumerate(["AFR","EUR","EAS","SAS","AMR"]):
        super_pops.extend([name] * N_PER_POP)
    super_arr = np.array(super_pops)
    n_hap = len(super_arr)

    # Assign each hap to a branch
    branch = np.zeros(n_hap, dtype=int)
    for s, name in enumerate(["AFR","EUR","EAS","SAS","AMR"]):
        idx = np.where(super_arr == name)[0]
        branch[idx] = rng.choice(N_BRANCHES, size=len(idx),
                                  p=BRANCH_FREQ[:, s])
    # SNV layout
    n_marker = N_BRANCHES * N_MARKER_PER_BRANCH
    n_backbone = N_BACKBONE
    n_noise = N_SNV - n_marker - n_backbone
    snv_branch = np.full(N_SNV, -1, dtype=int)   # -1 = backbone or noise
    for b in range(N_BRANCHES):
        snv_branch[b*N_MARKER_PER_BRANCH:(b+1)*N_MARKER_PER_BRANCH] = b
    # backbone markers: tag = N_BRANCHES (i.e., H4=3 already taken; assign tag=4 for backbone-spec)
    # Actually treat backbone differently: they have signal proportional to "carrier in H4"
    # but we also add a separate "backbone" group for general-purpose carriers
    bb_start = n_marker
    bb_end = bb_start + n_backbone
    # noise
    noise_start = bb_end

    # Build matrix
    M = np.zeros((n_hap, N_SNV), dtype=np.int8)
    # Branch markers
    for j in range(n_marker):
        b = snv_branch[j]
        carrier_p = np.where(branch == b, P_CARRIER_MATCH, P_CARRIER_MISMATCH)
        M[:, j] = (rng.random(n_hap) < carrier_p).astype(np.int8)
    # Backbone markers (any branch)
    for j in range(bb_start, bb_end):
        carrier_p = np.full(n_hap, P_BACKBONE_CARRIER)
        # boost slightly for H4
        carrier_p[branch == 3] = 0.6
        M[:, j] = (rng.random(n_hap) < carrier_p).astype(np.int8)
    # Noise: random MAF in [0.05, 0.45]
    for j in range(noise_start, N_SNV):
        maf = rng.uniform(0.05, 0.45)
        M[:, j] = (rng.random(n_hap) < maf).astype(np.int8)
    return M, super_arr, branch, snv_branch


def jaccard(a, b):
    if not a or not b: return 0.0
    return len(a & b) / len(a | b)


def run_pipeline(M, super_arr, branch, snv_branch, label="real"):
    """Run UMAP + Leiden + NMF k=4..6 + enrichment vs ground truth."""
    print(f"\n=== {label} ===")
    print(f"matrix: {M.shape}")
    Mfloat = M.astype(np.float32)

    # UMAP
    print("  UMAP...")
    reducer = umap.UMAP(n_neighbors=20, min_dist=0.05, metric="hamming",
                         random_state=0, n_components=2)
    xy = reducer.fit_transform(Mfloat)

    # Leiden
    print("  Leiden...")
    nn = NearestNeighbors(n_neighbors=20, metric="hamming"); nn.fit(Mfloat)
    knn = nn.kneighbors_graph(Mfloat, mode="connectivity")
    src, dst = knn.nonzero()
    edges = [(int(s), int(d)) for s, d in zip(src, dst) if s < d]
    g = ig.Graph(n=Mfloat.shape[0], edges=edges, directed=False)
    part = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=1.0, seed=0)
    leiden_labels = np.array(part.membership)
    n_leiden = len(set(leiden_labels))
    print(f"    {n_leiden} clusters")

    # NMF k=4 (matches ground truth)
    K_USE = 4
    print(f"  NMF k={K_USE}...")
    nmf = NMF(n_components=K_USE, init="nndsvdar", random_state=0,
               max_iter=400, tol=1e-4)
    W = nmf.fit_transform(Mfloat); H = nmf.components_

    # Per-component top 5% SNV signatures and ground-truth overlap
    n_snv = Mfloat.shape[1]
    top_n = max(1, int(NMF_TOP_FRAC * n_snv))
    overlap_mat = np.zeros((K_USE, N_BRANCHES))   # rows=NMF comp, cols=true branch
    for c in range(K_USE):
        top_idx = set(np.argsort(H[c])[-top_n:].tolist())
        for b in range(N_BRANCHES):
            true_set = set(np.where(snv_branch == b)[0].tolist())
            overlap_mat[c, b] = len(top_idx & true_set)
    print(f"  raw overlap matrix (NMF comp rows × true branch cols):")
    print(overlap_mat.astype(int))

    # Hungarian assignment: maximize total overlap (negate for minimization)
    cost = -overlap_mat
    row_ind, col_ind = linear_sum_assignment(cost)
    assignment = dict(zip(row_ind, col_ind))
    print(f"  Hungarian assignment NMF→branch: {assignment}")

    # Compute per-component stats (Jaccard, hypergeom p)
    rows = []
    for c in range(K_USE):
        top_idx = set(np.argsort(H[c])[-top_n:].tolist())
        for b in range(N_BRANCHES):
            true_set = set(np.where(snv_branch == b)[0].tolist())
            x = len(top_idx & true_set)
            K_succ = len(true_set)
            jac = jaccard(top_idx, true_set)
            p_val = float(hypergeom.sf(x - 1, n_snv, K_succ, top_n))
            rows.append({
                "label": label, "nmf_component": c, "true_branch": b,
                "overlap": x, "expected": top_n * K_succ / n_snv,
                "fold": x / (top_n * K_succ / n_snv) if K_succ > 0 else float("nan"),
                "jaccard": jac, "hyper_p": p_val,
                "assigned_to": assignment.get(c, -1),
            })

    # Branch composition by super-pop (real branch ground truth)
    bc = pd.DataFrame({"super_pop": super_arr, "true_branch": branch})
    bc_table = bc.groupby(["super_pop", "true_branch"]).size().unstack(fill_value=0)
    bc_freq = bc_table.div(bc_table.sum(axis=1), axis=0)
    print(f"  ground-truth branch frequencies by super-pop:")
    print(bc_freq.round(2).to_string())

    return {
        "label": label, "xy": xy, "leiden": leiden_labels,
        "W": W, "H": H, "overlap_mat": overlap_mat,
        "assignment": assignment, "rows": rows,
        "super_arr": super_arr, "branch": branch,
        "snv_branch": snv_branch, "bc_freq": bc_freq,
    }


def main():
    sys.stdout.reconfigure(line_buffering=True)
    t0 = time.time()

    # ---- Real ground-truth simulation ----
    M, super_arr, branch, snv_branch = simulate(rng_seed=42)
    real = run_pipeline(M, super_arr, branch, snv_branch, label="ground_truth")

    # ---- Negative control: shuffle population labels ----
    rng = np.random.default_rng(0)
    super_shuf = rng.permutation(super_arr.copy())
    # Keep the same hap matrix M (so SNV-branch correlations preserved internally),
    # but population labels are now random — UMAP/Leiden should still find branches
    # but they should NOT correlate with super_pops anymore. This tests that population-
    # private clustering reflects real population structure not artifact.
    neg = run_pipeline(M, super_shuf, branch, snv_branch, label="neg_ctrl_shuffled_pops")
    # In this case, branch composition by super_pop should be uniform.

    # ---- Pure-noise control: matrix without branch structure ----
    rng = np.random.default_rng(1)
    M_pure = (rng.random((M.shape[0], M.shape[1])) < 0.20).astype(np.int8)
    pure = run_pipeline(M_pure, super_arr, np.zeros(M.shape[0], dtype=int),
                         np.full(M.shape[1], -1, dtype=int), label="pure_noise")

    # ---- Combine results ----
    rows_all = real["rows"] + neg["rows"] + pure["rows"]
    df = pd.DataFrame(rows_all)
    df.to_csv(OUT / "sim_ground_truth_summary.csv", index=False)
    print(f"\nwrote {OUT/'sim_ground_truth_summary.csv'}")

    # Print summary tables
    for run in [real, neg, pure]:
        print(f"\n=== {run['label']}: NMF top-5% × true branch overlap (Jaccard) ===")
        sub = pd.DataFrame(run["rows"])
        piv = sub.pivot(index="nmf_component", columns="true_branch",
                          values="jaccard").round(3)
        print(piv.to_string())

    # ===== Figure: 3 rows × 4 cols =====
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":9,
        "savefig.bbox":"tight","pdf.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,
    })
    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1], hspace=0.35, wspace=0.3)

    super_palette = {"AFR":"#b04545","EUR":"#3a8c5a","EAS":"#5a4fcf",
                      "SAS":"#b5722d","AMR":"#7e3da7"}
    branch_palette = ["#d65454","#3070b8","#3a8c5a","#7e3da7"]

    for row, run in enumerate([real, neg, pure]):
        # A: UMAP colored by super-pop
        ax = fig.add_subplot(gs[row, 0])
        for s in ["AFR","EUR","EAS","SAS","AMR"]:
            m = run["super_arr"] == s
            ax.scatter(run["xy"][m, 0], run["xy"][m, 1], s=8,
                        color=super_palette[s], alpha=0.6,
                        edgecolor="none", label=s)
        ax.set_title(f"{run['label']} | A. UMAP super-pop", loc="left", fontsize=9)
        if row == 0:
            ax.legend(fontsize=7, loc="upper right")

        # B: UMAP colored by ground-truth branch
        ax = fig.add_subplot(gs[row, 1])
        for b in range(N_BRANCHES):
            m = run["branch"] == b
            if m.sum() == 0: continue
            ax.scatter(run["xy"][m, 0], run["xy"][m, 1], s=8,
                        color=branch_palette[b], alpha=0.6,
                        edgecolor="none", label=f"H{b+1}")
        ax.set_title(f"{run['label']} | B. UMAP true branch", loc="left", fontsize=9)
        if row == 0:
            ax.legend(fontsize=7, loc="upper right")

        # C: NMF top-5% × true branch overlap heatmap
        ax = fig.add_subplot(gs[row, 2])
        ovl = run["overlap_mat"]
        # Show as Jaccard
        jac_mat = np.zeros_like(ovl, dtype=float)
        n_snv = run["H"].shape[1]
        top_n = max(1, int(NMF_TOP_FRAC * n_snv))
        for c in range(ovl.shape[0]):
            for b in range(ovl.shape[1]):
                true_set_size = (run["snv_branch"] == b).sum()
                if true_set_size == 0:
                    jac_mat[c, b] = float("nan")
                else:
                    jac_mat[c, b] = ovl[c, b] / (top_n + true_set_size - ovl[c, b])
        im = ax.imshow(jac_mat, aspect="auto", cmap="Reds", vmin=0, vmax=1,
                        interpolation="nearest")
        ax.set_xticks(range(N_BRANCHES))
        ax.set_xticklabels([f"H{b+1}" for b in range(N_BRANCHES)])
        ax.set_yticks(range(ovl.shape[0]))
        ax.set_yticklabels([f"NMF c{c}" for c in range(ovl.shape[0])])
        for c in range(ovl.shape[0]):
            for b in range(ovl.shape[1]):
                v = jac_mat[c, b]
                if not np.isnan(v):
                    ax.text(b, c, f"{v:.2f}", ha="center", va="center",
                             fontsize=8, color="#fff" if v > 0.4 else "#111")
        fig.colorbar(im, ax=ax, fraction=0.04)
        ax.set_title(f"{run['label']} | C. NMF×true Jaccard", loc="left", fontsize=9)

        # D: branch composition by super-pop (ground-truth or shuffled)
        ax = fig.add_subplot(gs[row, 3])
        bc = run["bc_freq"]
        bottom = np.zeros(len(bc))
        x = np.arange(len(bc))
        for b in range(N_BRANCHES):
            if b not in bc.columns: continue
            ax.bar(x, bc[b].values, bottom=bottom, color=branch_palette[b],
                    edgecolor="#222", linewidth=0.4, label=f"H{b+1}")
            bottom += bc[b].values
        ax.set_xticks(x); ax.set_xticklabels(bc.index, fontsize=8)
        ax.set_ylabel("freq")
        ax.set_title(f"{run['label']} | D. branch by super-pop", loc="left", fontsize=9)
        if row == 0:
            ax.legend(fontsize=7, loc="upper right")
        ax.set_ylim(0, 1)

    fig.suptitle(
        "(f) Ground-truth simulation — NMF recovery of known latent branches",
        fontsize=11.5, y=1.005,
    )
    out = OUT / "figure_simulation_ground_truth"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\nwrote {out}")
    plt.close(fig)
    print(f"\nTotal: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
