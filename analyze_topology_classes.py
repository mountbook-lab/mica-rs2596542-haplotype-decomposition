"""Internal exploration: does Lange 2013's signed-correlation flip
decompose into a carrier-set topology class flip?

Pairs analysed:
  rs2596542 (anchor, MICA upstream)  ×  rs2244546 (HCP5 region, ~69 kb)
  rs2596542                          ×  rs9275572 (HLA-DQ region, ~1.3 Mb)

Lange 2013 Fig. 1A reported (signed Pearson r in HapMap3):
  rs2596542 × rs9275572: r = +0.27 (JPT), r = -0.12 (CEU)
  rs2596542 × rs2244546: r = +0.59 (JPT), r = -0.19 (CEU)

We re-render this in 1000G phase 3 (n=2,504, 26 sub-populations) and
add carrier-set topology class (DISJOINT / nested-A-in-B / nested-B-in-A
/ IDENTICAL / PARTIAL) plus the asymmetry index Δ = P(v|tag) − P(tag|v).

Output:
  results/per_pop_topology.csv
  figure_signed_r_vs_topology.{pdf,png}
"""

from __future__ import annotations

import os

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
VCF = REPO / "data/three_snps.vcf"
PANEL = Path(os.environ.get("MICA_DATA_DIR", "data")) / "integrated_call_samples_v3.20130502.ALL.panel"
OUT_DIR = REPO / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ANCHOR = ("rs2596542", 31_366_595)
PARTNERS = [
    ("rs2244546", 31_435_833,  "+0.59 JPT, -0.19 CEU"),
    ("rs9275572", 32_678_999,  "+0.27 JPT, -0.12 CEU"),
]

EAS_POPS = ["CHB", "JPT", "CHS", "CDX", "KHV"]
EUR_POPS = ["CEU", "TSI", "FIN", "GBR", "IBS"]


def parse_vcf(path: Path):
    """Return (samples, pos_to_alleles, pos_to_haplotypes).

    haplotypes is dict pos -> ndarray (n_haplotypes,) of 0/1 alleles
    in dbSNP forward-strand convention (REF=0, ALT=1).
    """
    samples = None
    haps_by_pos = {}
    alleles_by_pos = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
            # Keep only biallelic SNVs
            if alt.startswith("<") or "," in alt:
                continue
            gts = f[9:]
            n = len(gts)
            hap = np.empty(2 * n, dtype=np.int8)
            for i, g in enumerate(gts):
                # phased "a|b"
                a, b = g.split("|")
                hap[2 * i]     = int(a)
                hap[2 * i + 1] = int(b)
            haps_by_pos[pos] = hap
            alleles_by_pos[pos] = (ref, alt)
    return samples, alleles_by_pos, haps_by_pos


def topology_class(n11, n10, n01) -> str:
    """5-class operational definition (consistent with paper 1)."""
    if n11 == 0:
        return "DISJOINT"
    if n10 == 0 and n01 == 0:
        return "IDENTICAL"
    if n10 == 0 and n01 > 0:
        return "STRICT_A_IN_B"   # S_A subset S_B (POSITIVE-proxy)
    if n01 == 0 and n10 > 0:
        return "STRICT_B_IN_A"   # S_B subset S_A
    return "PARTIAL"


def analyse_pair(hap_a: np.ndarray, hap_b: np.ndarray, idx: np.ndarray):
    """Compute n11/n10/n01/n00, signed r, P(v|tag), P(tag|v), Δ, class
    on the haplotype subset selected by idx."""
    a = hap_a[idx]
    b = hap_b[idx]
    n = len(a)
    if n == 0:
        return None
    n11 = int(((a == 1) & (b == 1)).sum())
    n10 = int(((a == 1) & (b == 0)).sum())
    n01 = int(((a == 0) & (b == 1)).sum())
    n00 = int(((a == 0) & (b == 0)).sum())
    p_a = (n11 + n10) / n
    p_b = (n11 + n01) / n
    # Pearson r (signed). Standardise so REF=0 ALT=1 throughout.
    sa = a.std(); sb = b.std()
    if sa == 0 or sb == 0:
        r_signed = float("nan"); r_squared = float("nan")
    else:
        r_signed = float(((a - p_a) * (b - p_b)).mean() / (sa * sb))
        r_squared = r_signed ** 2
    p_b_given_a = (n11 / (n11 + n10)) if (n11 + n10) > 0 else float("nan")
    p_a_given_b = (n11 / (n11 + n01)) if (n11 + n01) > 0 else float("nan")
    delta = p_b_given_a - p_a_given_b
    cls = topology_class(n11, n10, n01)
    return dict(
        n=n, n11=n11, n10=n10, n01=n01, n00=n00,
        p_a=p_a, p_b=p_b,
        r_signed=r_signed, r_squared=r_squared,
        P_b_given_a=p_b_given_a, P_a_given_b=p_a_given_b,
        delta=delta, topology_class=cls,
    )


def main() -> None:
    sys.stdout.reconfigure(line_buffering=True)
    samples, alleles, haps = parse_vcf(VCF)
    print(f"loaded VCF: {len(samples)} samples")
    for pos, (ref, alt) in alleles.items():
        rs = next((rs for rs, p in [ANCHOR] + [(p[0], p[1]) for p in PARTNERS]
                    if p == pos), "?")
        print(f"  pos {pos:>10}  {ref}/{alt}   ({rs})")

    panel = pd.read_csv(PANEL, sep="\t")
    panel = panel[panel["sample"].isin(samples)]
    pop_of = dict(zip(panel["sample"], panel["pop"]))
    super_of = dict(zip(panel["sample"], panel["super_pop"]))

    # Build a haplotype-level pop array (2 * n_samples)
    hap_pop = np.empty(2 * len(samples), dtype=object)
    hap_super = np.empty(2 * len(samples), dtype=object)
    for i, s in enumerate(samples):
        hap_pop[2 * i]     = pop_of.get(s, "UNK")
        hap_pop[2 * i + 1] = pop_of.get(s, "UNK")
        hap_super[2 * i]     = super_of.get(s, "UNK")
        hap_super[2 * i + 1] = super_of.get(s, "UNK")

    # Anchor haplotypes
    anchor_pos = ANCHOR[1]
    if anchor_pos not in haps:
        raise SystemExit(f"anchor {ANCHOR[0]} ({anchor_pos}) not in VCF")
    anchor_hap = haps[anchor_pos]

    rows = []
    for pname, ppos, lange_note in PARTNERS:
        if ppos not in haps:
            print(f"  SKIP partner {pname} ({ppos}) -- not in VCF")
            continue
        partner_hap = haps[ppos]

        # Per-population
        for pop in EAS_POPS + EUR_POPS:
            idx = np.where(hap_pop == pop)[0]
            res = analyse_pair(anchor_hap, partner_hap, idx)
            if res is None: continue
            rows.append({"partner": pname, "pop": pop,
                          "super_pop": super_of[
                              [s for s in samples if pop_of.get(s) == pop][0]
                          ] if any(pop_of.get(s) == pop for s in samples) else "?",
                          "lange_note": lange_note,
                          **res})

        # Super-population pooled
        for sup in ["EAS", "EUR"]:
            idx = np.where(hap_super == sup)[0]
            res = analyse_pair(anchor_hap, partner_hap, idx)
            if res is None: continue
            rows.append({"partner": pname, "pop": f"{sup}_pooled",
                          "super_pop": sup,
                          "lange_note": lange_note,
                          **res})

    df = pd.DataFrame(rows)
    out_csv = OUT_DIR / "per_pop_topology.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}")

    # Print summary
    print("\n=== Summary ===")
    cols = ["partner", "pop", "n", "n11", "n10", "n01", "n00",
             "p_a", "p_b", "r_signed", "r_squared",
             "P_b_given_a", "P_a_given_b", "delta", "topology_class"]
    print(df[cols].to_string(
        index=False,
        formatters={
            "p_a": "{:.3f}".format,
            "p_b": "{:.3f}".format,
            "r_signed": "{:+.3f}".format,
            "r_squared": "{:.3f}".format,
            "P_b_given_a": "{:.3f}".format,
            "P_a_given_b": "{:.3f}".format,
            "delta": "{:+.3f}".format,
        },
    ))


if __name__ == "__main__":
    main()
