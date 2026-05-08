"""26-pop region-wide M/C/r/class scan around rs2596542.

For every biallelic SNV in chr6:31,116,595-31,616,595 (+/-250 kb),
compute per-population M, C, signed r, topology class across all 26
1000G phase 3 sub-populations + 5 super-pops (AFR, EUR, EAS, SAS, AMR).

Multi-super-pop C-flip definition:
  for each partner SNV, take median C within each super-pop (5 values)
  and count number of super-pop pairs with opposite-sign median C.
  "any_flip" = at least one pair flips. "n_flips" in [0, 10].

Outputs:
  results/region_per_pop_26.parquet
  results/region_per_partner_26.csv
"""

from __future__ import annotations

import gzip
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
VCF = REPO / "data/region_500kb.vcf.gz"
PANEL = Path("/home/yyamada1225/LVexome/demo_delta_mc/thousand_g/"
             "integrated_call_samples_v3.20130502.ALL.panel")
OUT = REPO / "results"

ANCHOR_POS = 31_366_595

ALL_26 = {
    "AFR": ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"],
    "EUR": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "EAS": ["CHB", "JPT", "CHS", "CDX", "KHV"],
    "SAS": ["GIH", "PJL", "BEB", "STU", "ITU"],
    "AMR": ["MXL", "PUR", "CLM", "PEL"],
}
KEEP_POPS = [p for super_ps in ALL_26.values() for p in super_ps]
SUPER_OF = {p: s for s, ps in ALL_26.items() for p in ps}
SUPER_POPS = list(ALL_26.keys())


def parse_vcf_streaming(path: Path):
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                yield ("SAMPLES", line.strip().split("\t")[9:]); continue
            f = line.rstrip("\n").split("\t")
            ref, alt = f[3], f[4]
            if (len(ref) != 1 or len(alt) != 1
                or alt.startswith("<") or "," in alt): continue
            pos = int(f[1])
            gts = f[9:]
            n = len(gts); hap = np.empty(2*n, dtype=np.int8)
            ok = True
            for i, g in enumerate(gts):
                if "|" not in g: ok = False; break
                a, b = g.split("|", 1)
                if a == "." or b == ".":
                    hap[2*i] = -1; hap[2*i+1] = -1
                else:
                    hap[2*i] = int(a); hap[2*i+1] = int(b)
            if not ok: continue
            yield (pos, ref, alt, hap)


def topology_class(n11, n10, n01) -> str:
    if n11 == 0: return "DISJOINT"
    if n10 == 0 and n01 == 0: return "IDENTICAL"
    if n10 == 0: return "STRICT_A_IN_B"
    if n01 == 0: return "STRICT_B_IN_A"
    return "PARTIAL"


def metrics(a: np.ndarray, b: np.ndarray):
    valid = (a >= 0) & (b >= 0)
    a = a[valid]; b = b[valid]
    n = len(a)
    if n == 0: return None
    n11 = int(((a == 1) & (b == 1)).sum())
    n10 = int(((a == 1) & (b == 0)).sum())
    n01 = int(((a == 0) & (b == 1)).sum())
    n00 = n - n11 - n10 - n01
    p_a = (n11 + n10) / n; p_b = (n11 + n01) / n
    if p_a == 0 or p_a == 1 or p_b == 0 or p_b == 1: return None
    P_AgB = n11 / (n11 + n01) if (n11 + n01) > 0 else float("nan")
    P_BgA = n11 / (n11 + n10) if (n11 + n10) > 0 else float("nan")
    Delta = P_AgB - P_BgA
    M = p_a - p_b
    C = Delta - M
    sa = a.std(); sb = b.std()
    r = float("nan") if (sa == 0 or sb == 0) else \
        float(((a - p_a) * (b - p_b)).mean() / (sa * sb))
    cls = topology_class(n11, n10, n01)
    return n, n11, p_a, p_b, M, C, Delta, r, cls


def main():
    sys.stdout.reconfigure(line_buffering=True)
    t0 = time.time()
    panel = pd.read_csv(PANEL, sep="\t")
    pop_of = dict(zip(panel["sample"], panel["pop"]))
    super_of = dict(zip(panel["sample"], panel["super_pop"]))

    print("loading VCF...")
    anchor_hap = None
    samples = None
    for rec in parse_vcf_streaming(VCF):
        if rec[0] == "SAMPLES":
            samples = rec[1]; continue
        pos = rec[0]
        if pos == ANCHOR_POS:
            anchor_hap = rec[3].copy()
            print(f"  found anchor at {pos}: {rec[1]}/{rec[2]}, "
                   f"alt freq global = "
                   f"{(anchor_hap[anchor_hap>=0]==1).mean():.3f}")
            break
    print(f"  loaded samples: {len(samples)}")

    hap_pop = np.empty(2*len(samples), dtype=object)
    hap_super = np.empty(2*len(samples), dtype=object)
    for i, s in enumerate(samples):
        hap_pop[2*i] = pop_of.get(s, "UNK"); hap_pop[2*i+1] = pop_of.get(s, "UNK")
        hap_super[2*i] = super_of.get(s, "UNK"); hap_super[2*i+1] = super_of.get(s, "UNK")

    # Per-pop carrier counts
    print("\nrs2596542-T (alt) carriers per pop:")
    for pop in KEEP_POPS:
        m = (hap_pop == pop) & (anchor_hap == 1)
        ntot = (hap_pop == pop).sum()
        print(f"  {pop:5s}: {m.sum():4d} / {ntot:4d} ({m.sum()/ntot*100:.1f}%)")

    pop_idx = {pop: np.where(hap_pop == pop)[0] for pop in KEEP_POPS}
    for sup in SUPER_POPS:
        pop_idx[f"{sup}_pooled"] = np.where(hap_super == sup)[0]

    print(f"\nscanning region for partners ({time.time()-t0:.0f}s)...")
    rows = []
    n_done = 0
    for rec in parse_vcf_streaming(VCF):
        if rec[0] == "SAMPLES": continue
        pos, ref, alt, hap = rec
        if pos == ANCHOR_POS: continue
        for pop, idx in pop_idx.items():
            res = metrics(anchor_hap[idx], hap[idx])
            if res is None: continue
            n, n11, p_a, p_b, M, C, Delta, r, cls = res
            rows.append((pos, pop, n, n11, p_a, p_b, M, C, Delta, r, cls))
        n_done += 1
        if n_done % 4000 == 0:
            print(f"  ...{n_done} partners ({time.time()-t0:.0f}s)")
    df = pd.DataFrame(rows, columns=[
        "partner_pos","pop","n","n11","p_a","p_b","M","C","Delta","r_signed","topology_class"])
    df.to_parquet(OUT / "region_per_pop_26.parquet")
    print(f"  wrote {OUT / 'region_per_pop_26.parquet'}  ({len(df)} rows)")

    # ---- per-partner summary across super-pops ----
    by = df[df["pop"].isin(KEEP_POPS)].copy()
    by["super"] = by["pop"].map(SUPER_OF)
    sup_med = by.groupby(["partner_pos","super"]).agg(
        median_C=("C","median"),
        median_M=("M","median"),
        median_r=("r_signed","median"),
        mode_class=("topology_class", lambda s: s.value_counts().idxmax()),
        any_disjoint=("topology_class", lambda s: (s == "DISJOINT").any()),
    ).reset_index()

    wide = sup_med.pivot(index="partner_pos", columns="super",
                          values=["median_C","median_M","median_r","mode_class","any_disjoint"]
                          ).reset_index()
    wide.columns = ["_".join(filter(None, [a, b])) for a, b in wide.columns]
    wide = wide.rename(columns={"partner_pos_": "partner_pos"})

    # multi-super-pop sign(C) flip
    cmat = wide[[f"median_C_{s}" for s in SUPER_POPS]].astype(float).values  # (n, 5)
    n_flips = np.zeros(len(wide), dtype=int)
    n_pairs = np.zeros(len(wide), dtype=int)
    for i, j in combinations(range(len(SUPER_POPS)), 2):
        ci = cmat[:, i]; cj = cmat[:, j]
        valid = ~(np.isnan(ci) | np.isnan(cj))
        n_pairs += valid.astype(int)
        flip = valid & ((np.sign(ci) * np.sign(cj)) < 0)
        n_flips += flip.astype(int)
    wide["n_super_pairs"] = n_pairs
    wide["n_super_C_flips"] = n_flips
    wide["any_C_flip"] = (n_flips >= 1)
    wide["all_pairs_flip"] = (n_flips == n_pairs) & (n_pairs > 0)

    # mode-class flip across super-pops (any pair differs)
    cls_cols = [f"mode_class_{s}" for s in SUPER_POPS]
    nset = wide[cls_cols].apply(lambda r: r.dropna().nunique(), axis=1)
    wide["n_super_class_distinct"] = nset
    wide["any_class_flip"] = (nset >= 2)

    wide.to_csv(OUT / "region_per_partner_26.csv", index=False)
    print(f"  wrote {OUT / 'region_per_partner_26.csv'}  ({len(wide)} partners)")

    # Summary
    n = len(wide)
    has5 = (wide["n_super_pairs"] == 10)
    print(f"\n=== Region summary (n={n} partners) ===")
    print(f"  partners with all 5 super-pops covered: {has5.sum()}")
    sub = wide[has5]
    print(f"  any C-flip among 10 super-pop pairs: "
          f"{sub['any_C_flip'].mean()*100:.1f}%")
    print(f"  >= 3 pairs flip:                     "
          f"{(sub['n_super_C_flips']>=3).mean()*100:.1f}%")
    print(f"  >= 5 pairs flip:                     "
          f"{(sub['n_super_C_flips']>=5).mean()*100:.1f}%")
    print(f"  any class flip across super-pops:    "
          f"{sub['any_class_flip'].mean()*100:.1f}%")

    # Distribution of n_flips
    print("\n  Distribution of n_super_C_flips (out of 10 possible):")
    nf_dist = sub["n_super_C_flips"].value_counts().sort_index()
    for k, v in nf_dist.items():
        print(f"    {int(k)} flips: {v:5d} ({v/has5.sum()*100:5.1f} %)")


if __name__ == "__main__":
    main()
