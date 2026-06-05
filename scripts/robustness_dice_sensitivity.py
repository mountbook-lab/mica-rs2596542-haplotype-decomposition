#!/usr/bin/env python3
"""Robustness Task 2b (optional, exploratory) — DICE immune-cell eQTL direction.

Follow-up to the eQTLGen direction-replication check: query DICE (Schmiedel 2018)
immune-cell cis-eQTL summary statistics, via the EBI eQTL Catalogue, for the c5
HLA-B / HLA-C candidate SNVs, and test whether the GTEx-defined c5 directions
(HLA-B up / HLA-C down) are directionally concordant across immune cell types.
This is EXPLORATORY sensitivity only, not primary evidence.

Access: eQTL Catalogue v2 (DICE = Schmiedel_2018, quant_method=ge, 15 cell types,
GRCh38). Per-dataset associations are read via the eQTL Catalogue REST API
(paginated over the small MHC window; remote tabix to the EBI FTP was attempted
first but returned intermittent libcurl errors in this environment). beta is the
ALT-allele effect (same convention as GTEx NES); variants are matched to
candidates by rsID and harmonised to the 1000G ALT allele.

Candidate source: results/robustness_eqtlgen_c5_direction.tsv (Task 2 output;
status == "harmonized" rows carry rsid + ref/alt + GTEx ALT-effect direction;
strand-ambiguous A/T and C/G SNVs were already dropped there). Task 1 clump
representatives are flagged where the optional --clump-tsv is given.

Standalone analysis; does NOT modify any existing figure / notebook / result.
If the catalogue or tabix is unreachable, it reports the failure point and stops.

Outputs:
  results/robustness_dice_direction_by_celltype.tsv
  results/robustness_dice_match_summary.tsv
"""
from __future__ import annotations
import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

API = "https://www.ebi.ac.uk/eqtl/api/v2"
STUDY_LABEL = "Schmiedel_2018"
GENE_ENSG = {"HLA-B": "ENSG00000234745", "HLA-C": "ENSG00000204525"}
ENSG_GENE = {v: k for k, v in GENE_ENSG.items()}
# GRCh38 window covering the c5 candidate SNVs (b37 ~31.14-31.32 Mb -> b38 +~32 kb)
WIN_START_B38, WIN_END_B38 = 31_100_000, 31_450_000
AMBIG = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def _get_json(url, tries=4):
    last = ""
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 400:            # past last page / bad offset: definitive
                raise
            last = f"HTTP {e.code}"
            time.sleep(2 + 2 * i)
        except Exception as e:           # transient network / rate-limit
            last = str(e)
            time.sleep(2 + 2 * i)
    raise RuntimeError(last)


def get_datasets():
    url = (f"{API}/datasets?study_label={STUDY_LABEL}&quant_method=ge&size=100")
    with urllib.request.urlopen(url, timeout=60) as r:
        data = json.load(r)
    return [(d["dataset_id"], d["sample_group"]) for d in data]


def fetch_region(dataset_id):
    """Return {(gene, rsid): (ref_b38, alt_b38, beta)} for HLA-B/HLA-C in the
    GRCh38 window, via the eQTL Catalogue associations API (paginated by
    position-ordered offset, early-stopping once past the window)."""
    out = {}
    unavailable = set()                      # gene not present in this dataset
    for gene, ensg in GENE_ENSG.items():
        start = 0
        for _ in range(40):                  # page cap (safety)
            url = (f"{API}/datasets/{dataset_id}/associations"
                   f"?molecular_trait_id={ensg}&size=1000&start={start}")
            try:
                data = _get_json(url)
            except Exception:
                if start == 0:               # gene absent from dataset (400/no results)
                    unavailable.add(gene)
                break                        # otherwise: past last page = end of data
            if not isinstance(data, list) or not data:
                break
            page_max = 0
            for x in data:
                pos = x.get("position", 0)
                page_max = max(page_max, pos)
                if pos < WIN_START_B38 or pos > WIN_END_B38:
                    continue
                rsid = x.get("rsid")
                ref = (x.get("ref") or "").upper()
                alt = (x.get("alt") or "").upper()
                beta = x.get("beta")
                if (rsid and str(rsid).startswith("rs") and beta is not None
                        and len(ref) == 1 and len(alt) == 1):
                    out[(gene, str(rsid))] = (ref, alt, float(beta))
            if page_max > WIN_END_B38:
                break
            start += 1000
            time.sleep(0.15)
    return out, unavailable


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--candidates",
                    default=None,
                    help="default <results-dir>/robustness_eqtlgen_c5_direction.tsv")
    ap.add_argument("--out-dir", default="results")
    args = ap.parse_args()
    RES = Path(args.results_dir)
    OUT = Path(args.out_dir)
    OUT.mkdir(parents=True, exist_ok=True)
    cand_path = Path(args.candidates) if args.candidates else \
        RES / "robustness_eqtlgen_c5_direction.tsv"

    if not cand_path.exists():
        sys.exit(f"[STOP] candidate file not found: {cand_path}\n"
                 "Run scripts/robustness_eqtlgen_sensitivity.py (Task 2) first.")

    cand = pd.read_csv(cand_path, sep="\t")
    cand = cand[(cand["status"] == "harmonized") & cand["rsid"].astype(str).str.startswith("rs")].copy()
    cand["gtex_dir"] = cand["gtex_dir"].astype(int)
    # candidate key: (gene, rsid) -> (ref, alt, gtex_dir)
    cand_map = {(r["gene"], r["rsid"]): (str(r["ref"]).upper(), str(r["alt"]).upper(),
                int(r["gtex_dir"])) for _, r in cand.iterrows()}
    n_cand = {g: int((cand["gene"] == g).sum()) for g in GENE_ENSG}
    print(f"candidates with rsid (harmonized): {n_cand}")

    # ---- fetch DICE datasets ----
    try:
        datasets = get_datasets()
    except Exception as e:
        sys.exit(f"[STOP] could not reach eQTL Catalogue API ({API}): {e}")
    print(f"DICE ge datasets: {len(datasets)}")

    rows = []
    for ds_id, cell in datasets:
        try:
            dice, unavail = fetch_region(ds_id)
        except Exception as e:
            print(f"  [skip] {cell} ({ds_id}): {e}")
            for gene in GENE_ENSG:
                rows.append({"cell_type": cell, "dataset_id": ds_id, "gene": gene,
                             "status": "fetch_failed", "n_matched": 0,
                             "n_harmonized": 0, "n_concordant": 0,
                             "concordance_rate": float("nan"), "binom_p_vs_0.5": float("nan")})
            continue
        for gene in GENE_ENSG:
            if gene in unavail:
                rows.append({"cell_type": cell, "dataset_id": ds_id, "gene": gene,
                             "status": "gene_not_in_dataset", "n_matched": 0,
                             "n_harmonized": 0, "n_concordant": 0,
                             "concordance_rate": float("nan"), "binom_p_vs_0.5": float("nan")})
                continue
            matched = harmon = conc = 0
            for (g, rsid), (cref, calt, gdir) in cand_map.items():
                if g != gene:
                    continue
                hit = dice.get((gene, rsid))
                if hit is None:
                    continue
                matched += 1
                dref, dalt, beta = hit
                if {dref, dalt} != {cref, calt}:
                    continue
                if (cref, calt) in AMBIG:
                    continue
                beta_alt = beta if dalt == calt else -beta
                if beta_alt == 0:
                    continue
                harmon += 1
                if np.sign(beta_alt) == gdir:
                    conc += 1
            rate = conc / harmon if harmon else float("nan")
            binom = stats.binomtest(conc, harmon, 0.5).pvalue if harmon else float("nan")
            rows.append({"cell_type": cell, "dataset_id": ds_id, "gene": gene,
                         "status": "ok", "n_matched": matched,
                         "n_harmonized": harmon, "n_concordant": conc,
                         "concordance_rate": round(rate, 3) if rate == rate else float("nan"),
                         "binom_p_vs_0.5": binom})
        print(f"  {cell}: HLA-B harm/conc, HLA-C harm/conc computed")

    det = pd.DataFrame(rows)
    det.to_csv(OUT / "robustness_dice_direction_by_celltype.tsv", sep="\t", index=False)

    # ---- summary per gene across cells (pooled harmonised SNV-cell observations) ----
    summ = []
    for gene in list(GENE_ENSG) + ["c5_overall"]:
        sub = det[det["status"] == "ok"]
        sub = sub if gene == "c5_overall" else sub[sub["gene"] == gene]
        H = int(sub["n_harmonized"].sum())
        C = int(sub["n_concordant"].sum())
        rate = C / H if H else float("nan")
        binom = stats.binomtest(C, H, 0.5).pvalue if H else float("nan")
        n_cells_ok = int((sub["n_harmonized"] > 0).sum())
        summ.append({"set": gene, "n_celltypes_with_data": n_cells_ok,
                     "n_harmonized_snv_cell_obs": H, "n_concordant": C,
                     "concordance_rate": round(rate, 3) if rate == rate else float("nan"),
                     "binom_p_vs_0.5": binom})
    summ = pd.DataFrame(summ)
    summ.to_csv(OUT / "robustness_dice_match_summary.tsv", sep="\t", index=False)

    print("\n=== DICE direction by cell type (ok rows) ===")
    print(det[det["status"] == "ok"].to_string(index=False))
    print("\n=== DICE summary per gene (pooled across cells) ===")
    print(summ.to_string(index=False))
    print("\nDICE = Schmiedel_2018 via eQTL Catalogue (GRCh38, quant=ge); beta = ALT effect, "
          "oriented to 1000G ALT; matched by rsID. EXPLORATORY sensitivity only.")


if __name__ == "__main__":
    main()
