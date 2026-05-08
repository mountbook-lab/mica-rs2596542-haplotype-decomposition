"""Within-branch carrier-set inclusion scan.

Hypothesis layer:
  rs2596542-T tags H1 union H2 union H3 (k=3 from cluster_anchor_haplotypes.py).
  For each branch Hk and each nearby SNV, classify the carrier-set
  relationship between (haplotype is in Hk) and (haplotype carries
  partner alt allele), within rs2596542-T carriers (universe = 675 hap).

Operational anchor for this scan: branch membership Hk.
Operational partner: every SNV in the +/-250 kb window with MAF >= 0.05
in carriers (already filtered: 7,303 SNVs).

Topology classes:
  IDENTICAL       -- S_Hk == S_partner  (perfect branch tag)
  STRICT_B_IN_A   -- S_partner subset S_Hk (sub-branch tag)
  STRICT_A_IN_B   -- S_Hk subset S_partner (Hk implies partner)
  DISJOINT        -- S_Hk and S_partner disjoint
  PARTIAL         -- otherwise

Output:
  results/within_branch_scan.csv         (one row per branch x partner_pos)
  results/within_branch_diagnostic.csv   (top tag SNVs per branch)
  results/within_branch_eas_vs_eur.csv   (EAS H_k vs EUR H_k allele-freq corr)
  results/figure_within_branch_scan.{pdf,png}
"""

from __future__ import annotations

import gzip
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
VCF = REPO / "data/region_500kb.vcf.gz"
PANEL = Path("/home/yyamada1225/LVexome/demo_delta_mc/thousand_g/"
             "integrated_call_samples_v3.20130502.ALL.panel")
BRANCH_PARQUET = REPO / "results/anchor_haplotype_branches.parquet"
OUT = REPO / "results"

ANCHOR_POS = 31_366_595
EAS_POPS = ["CHB", "JPT", "CHS", "CDX", "KHV"]
EUR_POPS = ["CEU", "TSI", "FIN", "GBR", "IBS"]
KEEP_POPS = EAS_POPS + EUR_POPS
MAF_FLOOR = 0.05
K = 3   # use k=3 branch labels


def load_haplotypes(path: Path):
    samples = None
    positions = []
    haps = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]; continue
            f = line.rstrip("\n").split("\t")
            ref, alt = f[3], f[4]
            if (len(ref) != 1 or len(alt) != 1
                or alt.startswith("<") or "," in alt): continue
            pos = int(f[1])
            gts = f[9:]
            n = len(gts)
            hap = np.empty(2*n, dtype=np.int8)
            ok = True
            for i, g in enumerate(gts):
                if "|" not in g: ok = False; break
                a, b = g.split("|", 1)
                if a == "." or b == ".":
                    hap[2*i] = -1; hap[2*i+1] = -1
                else:
                    hap[2*i] = int(a); hap[2*i+1] = int(b)
            if not ok: continue
            positions.append(pos); haps.append(hap)
    return samples, np.array(positions), np.array(haps, dtype=np.int8)


def topology_class(n11, n10, n01) -> str:
    if n11 == 0: return "DISJOINT"
    if n10 == 0 and n01 == 0: return "IDENTICAL"
    if n10 == 0: return "STRICT_A_IN_B"   # S_branch subset S_partner
    if n01 == 0: return "STRICT_B_IN_A"   # S_partner subset S_branch
    return "PARTIAL"


def main():
    sys.stdout.reconfigure(line_buffering=True)
    t0 = time.time()

    print("loading VCF...")
    samples, positions, haps = load_haplotypes(VCF)
    print(f"  {haps.shape[1]} haps, {len(positions)} SNVs ({time.time()-t0:.0f}s)")

    panel = pd.read_csv(PANEL, sep="\t")
    pop_of = dict(zip(panel["sample"], panel["pop"]))
    super_of = dict(zip(panel["sample"], panel["super_pop"]))
    hap_pop = np.empty(haps.shape[1], dtype=object)
    hap_super = np.empty(haps.shape[1], dtype=object)
    for i, s in enumerate(samples):
        hap_pop[2*i]   = pop_of.get(s, "UNK")
        hap_pop[2*i+1] = pop_of.get(s, "UNK")
        hap_super[2*i] = super_of.get(s, "UNK")
        hap_super[2*i+1] = super_of.get(s, "UNK")

    # rs2596542 carriers
    arow = np.where(positions == ANCHOR_POS)[0][0]
    is_carrier = (haps[arow] == 1) & np.isin(hap_pop, KEEP_POPS)
    carrier_idx = np.where(is_carrier)[0]
    print(f"  rs2596542-T carriers in 10 pops: {len(carrier_idx)}")

    # Branch labels (haplotype-level table from previous step)
    br_df = pd.read_parquet(BRANCH_PARQUET)
    # The previous script produced 675 haplotypes ordered the same way as
    # we encounter carriers here (it iterates the same VCF + same pop filter
    # in the same order). We'll re-derive ordering by pop sequence to align.
    # Sanity: lengths match.
    if len(br_df) != len(carrier_idx):
        raise SystemExit(
            f"branch table {len(br_df)} != carriers {len(carrier_idx)}"
        )
    branch_k = br_df[f"branch_k{K}"].values   # 1..K labels

    # Restrict matrices to carriers
    sub_haps = haps[:, carrier_idx]
    sub_pop = hap_pop[carrier_idx]
    sub_super = hap_super[carrier_idx]

    # Drop anchor row, drop SNVs with missing data, MAF filter within carriers
    mask = np.ones(len(positions), dtype=bool); mask[arow] = False
    sub_haps = sub_haps[mask]; sub_pos = positions[mask]
    valid = (sub_haps >= 0).all(axis=1)
    sub_haps = sub_haps[valid]; sub_pos = sub_pos[valid]
    af = sub_haps.mean(axis=1)
    info = (af >= MAF_FLOOR) & (af <= 1 - MAF_FLOOR)
    sub_haps = sub_haps[info]; sub_pos = sub_pos[info]
    print(f"  informative SNVs: {len(sub_pos)}")

    # ===== Within-branch classification =====
    print(f"\nclassifying branch x partner (k={K})...")
    rows = []
    for k in range(1, K + 1):
        in_branch = (branch_k == k)
        n_branch = int(in_branch.sum())
        for j, pos in enumerate(sub_pos):
            partner = sub_haps[j]
            n11 = int((in_branch & (partner == 1)).sum())
            n10 = int((in_branch & (partner == 0)).sum())
            n01 = int((~in_branch & (partner == 1)).sum())
            n00 = int((~in_branch & (partner == 0)).sum())
            cls = topology_class(n11, n10, n01)
            rows.append({
                "branch": f"H{k}", "partner_pos": int(pos),
                "n_branch": n_branch, "n11": n11, "n10": n10,
                "n01": n01, "n00": n00,
                "p_branch": n_branch / (n_branch + n00 + n01),
                "p_partner": (n11 + n01) / (n11 + n10 + n01 + n00),
                "P_partner_given_branch": n11/(n11+n10) if (n11+n10) > 0 else float("nan"),
                "P_branch_given_partner": n11/(n11+n01) if (n11+n01) > 0 else float("nan"),
                "topology_class": cls,
            })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "within_branch_scan.csv", index=False)
    print(f"  wrote {OUT/'within_branch_scan.csv'} ({len(df)} rows)")

    # ===== Per-branch class distribution =====
    print()
    pivot = df.groupby(["branch", "topology_class"]).size().unstack(fill_value=0)
    print(pivot.to_string())
    print()
    pivot_pct = (pivot.div(pivot.sum(axis=1), axis=0) * 100).round(1)
    print("(% of partners per class)")
    print(pivot_pct.to_string())

    # ===== Diagnostic SNVs per branch =====
    print("\n=== Diagnostic SNVs (IDENTICAL = perfect branch tag) ===")
    diag = []
    for k in range(1, K + 1):
        sub = df[(df["branch"] == f"H{k}") &
                 (df["topology_class"] == "IDENTICAL")]
        print(f"  H{k}: {len(sub)} IDENTICAL partners")
        if len(sub):
            print(sub.head(5)[["partner_pos","n11","n_branch"]].to_string(index=False))
        sub2 = df[(df["branch"] == f"H{k}") &
                  (df["topology_class"] == "STRICT_B_IN_A")]
        print(f"      + {len(sub2)} STRICT_B_IN_A (S_partner subset S_Hk)")
        sub3 = df[(df["branch"] == f"H{k}") &
                  (df["topology_class"] == "STRICT_A_IN_B")]
        print(f"      + {len(sub3)} STRICT_A_IN_B (S_Hk subset S_partner)")
        diag.append(pd.concat([sub, sub2.head(20)]).head(50).assign(branch=f"H{k}"))
    diag_df = pd.concat(diag, ignore_index=True) if diag else pd.DataFrame()
    if len(diag_df):
        diag_df.to_csv(OUT / "within_branch_diagnostic.csv", index=False)
        print(f"\n  wrote {OUT/'within_branch_diagnostic.csv'}")

    # ===== EAS vs EUR within-branch agreement =====
    print("\n=== EAS vs EUR within-branch allele-frequency agreement ===")
    eas_in = np.isin(sub_super, ["EAS"])
    eur_in = np.isin(sub_super, ["EUR"])
    ev_rows = []
    for k in range(1, K + 1):
        in_branch = (branch_k == k)
        m_eas = in_branch & eas_in
        m_eur = in_branch & eur_in
        n_eas = m_eas.sum(); n_eur = m_eur.sum()
        if n_eas < 5 or n_eur < 5:
            continue
        af_eas = sub_haps[:, m_eas].mean(axis=1)
        af_eur = sub_haps[:, m_eur].mean(axis=1)
        # restrict to SNVs informative in either pop
        keep = ((af_eas > 0.02) | (af_eur > 0.02)) & ((af_eas < 0.98) | (af_eur < 0.98))
        if keep.sum() < 50:
            continue
        a, b = af_eas[keep], af_eur[keep]
        if a.std() == 0 or b.std() == 0:
            corr = float("nan")
        else:
            corr = float(np.corrcoef(a, b)[0, 1])
        ev_rows.append({
            "branch": f"H{k}", "n_EAS": int(n_eas), "n_EUR": int(n_eur),
            "n_snv_kept": int(keep.sum()),
            "pearson_AF_EAS_vs_EUR": corr,
            "median_abs_diff": float(np.median(np.abs(a - b))),
        })
        print(f"  H{k}: n_EAS={n_eas} n_EUR={n_eur} "
               f"pearson(AF_EAS, AF_EUR) = {corr:.3f}  "
               f"median |diff| = {np.median(np.abs(a - b)):.3f}")
    pd.DataFrame(ev_rows).to_csv(OUT / "within_branch_eas_vs_eur.csv", index=False)

    # ===== Figure: spatial scan of branch class =====
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.bbox": "tight", "pdf.fonttype": 42,
    })
    # Landmarks (GRCh37) for context
    LANDMARKS = [
        (31_366_595, "rs2596542\n(anchor, MICA up.)", "#222"),
        (31_435_833, "rs2244546\n(HCP5 region)",      "#444"),
        (31_475_084, "rs11509487\n(MICB)",            "#b04545"),
    ]
    fig, axes = plt.subplots(K, 1, figsize=(13, 7.5), sharex=True)
    cls_color = {
        "IDENTICAL": "#7e3da7",         # purple — perfect tag
        "STRICT_B_IN_A": "#3a8c5a",     # green — partner subset of branch
        "STRICT_A_IN_B": "#b5722d",     # brown — branch subset of partner
        "DISJOINT": "#4078a3",          # blue
        "PARTIAL": "#c8c8c8",           # grey
    }
    for ax, k in zip(axes, range(1, K + 1)):
        sub = df[df["branch"] == f"H{k}"]
        for cls, col in cls_color.items():
            m = sub["topology_class"] == cls
            ax.scatter(sub.loc[m, "partner_pos"]/1e6,
                        np.full(m.sum(), {
                            "IDENTICAL": 4, "STRICT_B_IN_A": 3,
                            "STRICT_A_IN_B": 2, "DISJOINT": 1,
                            "PARTIAL": 0,
                        }[cls]),
                        s=8 if cls != "PARTIAL" else 3,
                        color=col, alpha=0.7 if cls != "PARTIAL" else 0.3,
                        edgecolor="none", label=cls)
        for lpos, lname, lcolor in LANDMARKS:
            ax.axvline(lpos/1e6, color=lcolor, lw=0.7, ls="--", alpha=0.85)
            if k == 1:
                ax.text(lpos/1e6, 4.7, lname, fontsize=6.5,
                         color=lcolor, ha="center", va="bottom")
        ax.set_yticks([0, 1, 2, 3, 4])
        ax.set_yticklabels(["PARTIAL", "DISJOINT",
                             "STRICT_A_IN_B\n(branch ⊂ partner)",
                             "STRICT_B_IN_A\n(partner ⊂ branch)",
                             "IDENTICAL"], fontsize=7)
        ax.set_title(f"H{k} (n={int((branch_k==k).sum())} carriers)",
                      loc="left", fontsize=9.5)
        ax.set_ylim(-0.5, 5.2 if k == 1 else 4.5)
        ax.grid(axis="x", linestyle=":", alpha=0.3)
        if ax is axes[0]:
            ax.legend(fontsize=7, loc="upper right", ncol=5)
    axes[-1].set_xlabel("chromosomal position (Mb, GRCh37)")
    fig.suptitle(
        "Within-branch carrier-set inclusion scan — rs2596542-T carriers",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    out = OUT / "figure_within_branch_scan"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"\n  wrote {out}.png and .pdf")
    plt.close(fig)

    print(f"\nTotal: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
