"""Region-wide M/C/r/class scan around rs2596542.

For every biallelic SNV in chr6:31,116,595-31,616,595 (±250 kb around
rs2596542), compute per-population M, C, signed r, topology class
in 5 EAS pops (CHB,JPT,CHS,CDX,KHV) and 5 EUR pops (CEU,TSI,FIN,GBR,IBS),
plus EAS_pooled and EUR_pooled.

Output:
  results/region_per_pop.parquet  (compact long format)
  results/region_per_partner.csv  (per-partner summary across pops)
"""

from __future__ import annotations

import os

import gzip
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
VCF = REPO / "data/region_500kb.vcf.gz"
PANEL = Path(os.environ.get("MICA_DATA_DIR", "data")) / "integrated_call_samples_v3.20130502.ALL.panel"
OUT = REPO / "results"
OUT.mkdir(exist_ok=True)

ANCHOR_POS = 31_366_595   # rs2596542
EAS_POPS = ["CHB", "JPT", "CHS", "CDX", "KHV"]
EUR_POPS = ["CEU", "TSI", "FIN", "GBR", "IBS"]
KEEP_POPS = EAS_POPS + EUR_POPS


def parse_vcf_streaming(path: Path):
    """Yield (pos, ref, alt, hap_array) for each biallelic SNV; first
    yield is ('SAMPLES', samples)."""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                yield ("SAMPLES", samples)
                continue
            f = line.rstrip("\n").split("\t")
            ref, alt = f[3], f[4]
            if len(ref) != 1 or len(alt) != 1:
                continue
            if alt.startswith("<") or "," in alt:
                continue
            pos = int(f[1])
            gts = f[9:]
            n = len(gts)
            hap = np.empty(2 * n, dtype=np.int8)
            ok = True
            for i, g in enumerate(gts):
                # phased "a|b" only
                if "|" not in g:
                    ok = False
                    break
                a, b = g.split("|", 1)
                if a == "." or b == ".":
                    hap[2 * i] = -1
                    hap[2 * i + 1] = -1
                else:
                    hap[2 * i] = int(a)
                    hap[2 * i + 1] = int(b)
            if not ok:
                continue
            yield (pos, ref, alt, hap)


def topology_class(n11, n10, n01) -> str:
    if n11 == 0:
        return "DISJOINT"
    if n10 == 0 and n01 == 0:
        return "IDENTICAL"
    if n10 == 0:
        return "STRICT_A_IN_B"
    if n01 == 0:
        return "STRICT_B_IN_A"
    return "PARTIAL"


def metrics(a: np.ndarray, b: np.ndarray):
    valid = (a >= 0) & (b >= 0)
    a = a[valid]; b = b[valid]
    n = len(a)
    if n == 0:
        return None
    n11 = int(((a == 1) & (b == 1)).sum())
    n10 = int(((a == 1) & (b == 0)).sum())
    n01 = int(((a == 0) & (b == 1)).sum())
    n00 = n - n11 - n10 - n01
    p_a = (n11 + n10) / n
    p_b = (n11 + n01) / n
    if p_a == 0 or p_a == 1 or p_b == 0 or p_b == 1:
        return None
    P_AgB = n11 / (n11 + n01) if (n11 + n01) > 0 else float("nan")
    P_BgA = n11 / (n11 + n10) if (n11 + n10) > 0 else float("nan")
    Delta = P_AgB - P_BgA
    M = p_a - p_b
    C = Delta - M
    sa = a.std(); sb = b.std()
    if sa == 0 or sb == 0:
        r_signed = float("nan")
    else:
        r_signed = float(((a - p_a) * (b - p_b)).mean() / (sa * sb))
    cls = topology_class(n11, n10, n01)
    return n, n11, p_a, p_b, M, C, Delta, r_signed, cls


def main():
    sys.stdout.reconfigure(line_buffering=True)
    t0 = time.time()
    panel = pd.read_csv(PANEL, sep="\t")
    pop_of = dict(zip(panel["sample"], panel["pop"]))
    super_of = dict(zip(panel["sample"], panel["super_pop"]))

    # First pass: find anchor + collect samples
    anchor_hap = None
    samples = None
    print("scanning VCF...")
    pos_list = []
    n_snv = 0
    for rec in parse_vcf_streaming(VCF):
        if rec[0] == "SAMPLES":
            samples = rec[1]
            continue
        pos = rec[0]
        if pos == ANCHOR_POS:
            anchor_hap = rec[3].copy()
            print(f"  found anchor at {pos}: ref={rec[1]} alt={rec[2]}, "
                  f"alt freq = {(anchor_hap[anchor_hap>=0]==1).mean():.3f}")
        n_snv += 1
        pos_list.append(pos)
    if anchor_hap is None:
        raise SystemExit(f"anchor at {ANCHOR_POS} not found in VCF")
    print(f"  {n_snv} biallelic SNVs in window")

    # Build hap-level pop arrays
    hap_pop = np.empty(2 * len(samples), dtype=object)
    hap_super = np.empty(2 * len(samples), dtype=object)
    for i, s in enumerate(samples):
        hap_pop[2*i]     = pop_of.get(s, "UNK")
        hap_pop[2*i+1]   = pop_of.get(s, "UNK")
        hap_super[2*i]   = super_of.get(s, "UNK")
        hap_super[2*i+1] = super_of.get(s, "UNK")

    # Pre-compute index sets
    pop_idx = {pop: np.where(hap_pop == pop)[0] for pop in KEEP_POPS}
    pop_idx["EAS_pooled"] = np.where(hap_super == "EAS")[0]
    pop_idx["EUR_pooled"] = np.where(hap_super == "EUR")[0]

    # Second pass: compute metrics for every (partner × pop)
    rows = []
    n_done = 0
    for rec in parse_vcf_streaming(VCF):
        if rec[0] == "SAMPLES":
            continue
        pos, ref, alt, hap = rec
        if pos == ANCHOR_POS:
            continue
        for pop, idx in pop_idx.items():
            res = metrics(anchor_hap[idx], hap[idx])
            if res is None: continue
            n, n11, p_a, p_b, M, C, Delta, r_signed, cls = res
            rows.append((pos, pop, n, n11, p_a, p_b, M, C, Delta, r_signed, cls))
        n_done += 1
        if n_done % 2000 == 0:
            print(f"  ...{n_done} partners ({time.time()-t0:.0f}s)")

    print(f"  {n_done} partner SNVs processed in {time.time()-t0:.0f}s")
    df = pd.DataFrame(rows, columns=[
        "partner_pos","pop","n","n11","p_a","p_b","M","C","Delta","r_signed","topology_class"])
    df.to_parquet(OUT / "region_per_pop.parquet")
    print(f"  wrote {OUT / 'region_per_pop.parquet'}  ({len(df)} rows)")

    # Per-partner summary: sign(C) by EAS vs EUR sub-pops, class flip detection
    by = df[df["pop"].isin(KEEP_POPS)].copy()
    by["super_pop"] = by["pop"].apply(lambda p: "EAS" if p in EAS_POPS else "EUR")
    grp = by.groupby(["partner_pos","super_pop"]).agg(
        n_pos_C=("C", lambda s: (s > 0).sum()),
        n_neg_C=("C", lambda s: (s < 0).sum()),
        median_C=("C", "median"),
        median_M=("M", "median"),
        median_r=("r_signed", "median"),
        classes=("topology_class", lambda s: ",".join(sorted(set(s)))),
    ).reset_index()

    # Pivot to wide
    wide = grp.pivot(index="partner_pos", columns="super_pop", values=[
        "n_pos_C","n_neg_C","median_C","median_M","median_r","classes"]).reset_index()
    wide.columns = ["_".join(filter(None, [a, b])) for a, b in wide.columns]
    wide = wide.rename(columns={"partner_pos_": "partner_pos"})
    wide.to_csv(OUT / "region_per_partner.csv", index=False)
    print(f"  wrote {OUT / 'region_per_partner.csv'}  ({len(wide)} partners)")

    # Quick stats
    n_total = len(wide)
    has_eas = wide["median_C_EAS"].notna()
    has_eur = wide["median_C_EUR"].notna()
    both = has_eas & has_eur
    n_both = both.sum()
    eas_pos = (wide.loc[both, "median_C_EAS"] > 0)
    eas_neg = (wide.loc[both, "median_C_EAS"] < 0)
    eur_pos = (wide.loc[both, "median_C_EUR"] > 0)
    eur_neg = (wide.loc[both, "median_C_EUR"] < 0)
    flip_pos_neg = ((eas_pos) & (eur_neg)).sum()
    flip_neg_pos = ((eas_neg) & (eur_pos)).sum()
    same_pos = ((eas_pos) & (eur_pos)).sum()
    same_neg = ((eas_neg) & (eur_neg)).sum()

    # Class-flip detection: any partner where the set of classes differs
    # between EAS and EUR
    class_flip = (wide["classes_EAS"] != wide["classes_EUR"]).sum()

    print()
    print(f"=== Region summary ({n_both} partners with data in both EAS+EUR) ===")
    print(f"  median(C) sign pattern:")
    print(f"    EAS+ / EUR+ (no flip): {same_pos:6d} ({same_pos/n_both*100:.1f} %)")
    print(f"    EAS- / EUR- (no flip): {same_neg:6d} ({same_neg/n_both*100:.1f} %)")
    print(f"    EAS+ / EUR- (flip):    {flip_pos_neg:6d} ({flip_pos_neg/n_both*100:.1f} %)")
    print(f"    EAS- / EUR+ (flip):    {flip_neg_pos:6d} ({flip_neg_pos/n_both*100:.1f} %)")
    print(f"  any class change EAS vs EUR (incl. order): {class_flip} ({class_flip/n_total*100:.1f} %)")


if __name__ == "__main__":
    main()
