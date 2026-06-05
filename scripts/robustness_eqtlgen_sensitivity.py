#!/usr/bin/env python3
"""Robustness Task 2 — independent-eQTL direction replication for the c5 axis.

GTEx HLA expression can be affected by reference-mapping bias. As a best-effort
sensitivity check (NOT an HLA-aware re-quantification), we test whether the c5
signature direction (HLA-B up / HLA-C down, defined by GTEx NES sign) reproduces
in an independent dataset and pipeline: eQTLGen whole-blood cis-eQTL meta-analysis
(GRCh37). We compare the sign of the ALT-allele effect between GTEx (NES) and
eQTLGen (Z-score), per SNV, after allele harmonisation.

This is a STANDALONE robustness analysis. It does NOT modify or regenerate any
existing figure, table, notebook, or result file.

eQTLGen extract (prepared by streaming the FDR<0.05 significant cis-eQTL file and
keeping GeneSymbol in {HLA-B, HLA-C}; columns: Pvalue SNP SNPChr SNPPos
AssessedAllele OtherAllele Zscore Gene GeneSymbol GeneChr GenePos ...):
  source: https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/
          2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
  Z-score sign is the effect of AssessedAllele on expression.

Inputs (GRCh37; argparse overridable):
  results/_axis_nes_full.csv             c5 signature x HLA-B/HLA-C GTEx NES (b37)
  data/region_500kb.vcf.gz               REF/ALT alleles per b37 position
  results/eqtlgen_HLA-B_HLA-C_cis.tsv    eQTLGen extract (see above)

Outputs:
  results/robustness_eqtlgen_c5_direction.tsv
  results/robustness_eqtlgen_match_summary.tsv
"""
from __future__ import annotations
import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

NES_TISSUE = "Whole_Blood"
AMBIG = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def load_vcf_alleles(vcf_path, want_pos):
    """Return {b37_pos: (REF, ALT)} for biallelic SNVs at the wanted positions."""
    alleles = {}
    want = set(want_pos)
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 5)
            pos = int(f[1])
            if pos not in want:
                continue
            ref, alt = f[3], f[4]
            if len(ref) == 1 and len(alt) == 1 and alt not in (".", "<"):
                alleles[pos] = (ref.upper(), alt.upper())
    return alleles


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--vcf", default="data/region_500kb.vcf.gz")
    ap.add_argument("--eqtlgen", default=None,
                    help="eQTLGen HLA-B/HLA-C extract tsv "
                         "(default: <results-dir>/eqtlgen_HLA-B_HLA-C_cis.tsv)")
    ap.add_argument("--out-dir", default="results")
    args = ap.parse_args()
    RES = Path(args.results_dir)
    OUT = Path(args.out_dir)
    OUT.mkdir(parents=True, exist_ok=True)
    eqtlgen_path = Path(args.eqtlgen) if args.eqtlgen else RES / "eqtlgen_HLA-B_HLA-C_cis.tsv"

    # ---- c5 GTEx signature directions (HLA-B / HLA-C) ----
    nes = pd.read_csv(RES / "_axis_nes_full.csv")
    c5 = nes[(nes["component"] == "c5") & (nes["tissue"] == NES_TISSUE)
             & (nes["gene"].isin(["HLA-B", "HLA-C"]))].copy()
    c5["b37_pos"] = c5["b37_pos"].astype(int)
    # one GTEx direction per (gene, pos): use the listed NES (dedup if needed)
    c5 = c5.groupby(["gene", "b37_pos"], as_index=False)["NES"].mean()

    # ---- REF/ALT from the 1000G VCF (b37) ----
    alleles = load_vcf_alleles(args.vcf, c5["b37_pos"].tolist())

    # ---- eQTLGen extract ----
    eg = pd.read_csv(eqtlgen_path, sep="\t")
    eg = eg[eg["SNPChr"] == 6]
    eg_map = {}   # (GeneSymbol, SNPPos) -> (rsid, assessed, other, Z)
    for _, r in eg.iterrows():
        eg_map[(r["GeneSymbol"], int(r["SNPPos"]))] = (
            r["SNP"], str(r["AssessedAllele"]).upper(),
            str(r["OtherAllele"]).upper(), float(r["Zscore"]))

    rows = []
    for _, r in c5.iterrows():
        gene, pos, gtex_nes = r["gene"], int(r["b37_pos"]), float(r["NES"])
        rec = {"gene": gene, "b37_pos": pos, "gtex_nes": round(gtex_nes, 4),
               "gtex_dir": int(np.sign(gtex_nes)), "rsid": "", "ref": "", "alt": "",
               "eqtlgen_assessed": "", "eqtlgen_other": "", "eqtlgen_Z": np.nan,
               "eqtlgen_Z_alt_oriented": np.nan, "concordant": np.nan,
               "status": ""}
        if pos not in alleles:
            rec["status"] = "drop_no_biallelic_ref_alt"
            rows.append(rec); continue
        ref, alt = alleles[pos]
        rec["ref"], rec["alt"] = ref, alt
        if (ref, alt) in AMBIG:
            rec["status"] = "drop_ambiguous_AT_CG"
            rows.append(rec); continue
        if (gene, pos) not in eg_map:
            rec["status"] = "not_in_eqtlgen_significant"
            rows.append(rec); continue
        rsid, assessed, other, z = eg_map[(gene, pos)]
        rec.update({"rsid": rsid, "eqtlgen_assessed": assessed,
                    "eqtlgen_other": other, "eqtlgen_Z": round(z, 4)})
        if {assessed, other} != {ref, alt}:
            rec["status"] = "drop_allele_mismatch"
            rows.append(rec); continue
        # orient eQTLGen Z to the ALT allele (GTEx NES is ALT-allele effect)
        z_alt = z if assessed == alt else -z
        rec["eqtlgen_Z_alt_oriented"] = round(z_alt, 4)
        rec["concordant"] = bool(np.sign(z_alt) == np.sign(gtex_nes))
        rec["status"] = "harmonized"
        rows.append(rec)

    det = pd.DataFrame(rows)
    det.to_csv(OUT / "robustness_eqtlgen_c5_direction.tsv", sep="\t", index=False)

    # ---- summary per gene + overall ----
    summ = []
    for gene in ["HLA-B", "HLA-C", "c5_overall"]:
        sub = det if gene == "c5_overall" else det[det["gene"] == gene]
        n_cand = len(sub)
        n_matched = int((sub["status"] != "not_in_eqtlgen_significant").sum()
                        - sub["status"].isin(["drop_no_biallelic_ref_alt",
                                              "drop_ambiguous_AT_CG"]).sum())
        harm = sub[sub["status"] == "harmonized"]
        n_harm = len(harm)
        n_conc = int(harm["concordant"].sum())
        rate = n_conc / n_harm if n_harm else float("nan")
        binom = stats.binomtest(n_conc, n_harm, 0.5).pvalue if n_harm else float("nan")
        drops = sub["status"].value_counts().to_dict()
        summ.append({
            "set": gene,
            "n_candidate_snv": n_cand,
            "n_harmonized": n_harm,
            "n_concordant": n_conc,
            "concordance_rate": round(rate, 3) if rate == rate else float("nan"),
            "binom_p_vs_0.5": binom,
            "n_drop_no_alleles": drops.get("drop_no_biallelic_ref_alt", 0),
            "n_drop_ambiguous": drops.get("drop_ambiguous_AT_CG", 0),
            "n_not_in_eqtlgen": drops.get("not_in_eqtlgen_significant", 0),
            "n_drop_allele_mismatch": drops.get("drop_allele_mismatch", 0),
        })
    summ = pd.DataFrame(summ)
    summ.to_csv(OUT / "robustness_eqtlgen_match_summary.tsv", sep="\t", index=False)

    print("=== eQTLGen c5 direction replication (GTEx NES vs eQTLGen Z, ALT-oriented) ===")
    print(summ.to_string(index=False))
    print("\nstatus breakdown:")
    print(det["status"].value_counts().to_string())
    print(f"\nNES tissue = {NES_TISSUE}. eQTLGen = whole-blood cis-eQTL meta-analysis "
          "(FDR<0.05, GRCh37). Z-score sign = AssessedAllele effect, oriented to ALT.")
    print("This is independent-dataset DIRECTION replication, NOT HLA-aware "
          "re-quantification (HLApers/arcasHLA gold standard not performed).")


if __name__ == "__main__":
    main()
