# Targeted add-on audit — revised Figure 2, Supplementary Figure S7, Supplementary Table S2, and numbering

Date: 2026-05-15
Auditor: Claude Code (hybrid pass; not a full re-audit)
Repo root: `the working repository/`
Working dir: `mica_flipflop/`

This is an **add-on** audit limited to the components that were added
or revised after the prior full audit. The five scope items are:

1. Revised Figure 2 — k = 8 NMF component composition among rs2596542-T carrier haplotypes.
2. Supplementary Figure S7 — COMT Val158Met (rs4680) reference-panel replication and NMF extension.
3. Supplementary Table S2 — top-signature SNVs defining the NMF components.
4. Supplementary figure / table numbering and caption consistency.
5. GitHub-ready update summary (see `reports/github_update_summary_final.md`).

The scientific narrative to preserve is unchanged:

- rs2596542-T carrier haplotypes are decomposed by NMF at k = 8.
- Axis I = c4 / c6, MICA cis-regulatory, population-conserved.
- Axis II = c5, HLA-B / HLA-C regulatory, population-variable, signed-LD reversal enriched.
- c5 shows zero MICA eQTL overlap across tested thresholds.
- LIRI-JP validates a tumor-retained HLA-C downregulation arm.
- COMT Val158Met is only an external proof-of-concept (Supplementary S7), not a second main locus.
- NMF components are latent haplotypic SNV signatures, not ancestry components.

---

## 1. Files inspected

| File | Role |
|---|---|
| `figures/main/figure_2_nmf_k8_population_composition.{pdf,png}` | Revised main Fig. 2 |
| `figures/supplementary/figure_2_nmf_k8_26pop.{pdf,png}` | 26-pop companion |
| `figures/supplementary/figure_S4_comt_flipflop.{pdf,png}` | New Supp. Fig. S7 (COMT) |
| `figures/supplementary/figure_S_coarse_carrier_space_k2_k5.{pdf,png}` | Archived legacy Fig 2 (k=2..5) |
| `figure_2_nmf_k8_population_composition.py` | Source script for revised Fig. 2 |
| `figure_S4_comt_flipflop.py` | Source script for Supp. Fig. S7 |
| `supplementary/tables/table_S2_component_top_snvs.csv` | Top-signature SNV table |
| `supplementary/SUPPLEMENTARY_INDEX.md` | Pre-AJHG supplementary index |
| `tables/figure_2_component_population_summary.csv` | Fig. 2 companion CSV |
| `docs/FIGURE_NUMBERING.md` | Locked figure-numbering doc |
| `docs/manuscript_revisions_v2.md` | v2 manuscript revisions |
| `reports/figure_2_revision_audit.md` | Prior revision audit (Fig. 2) |
| `README.md` (`the working repository/README.md`) | Repo README |
| `supplementary/figures/` | Pre-AJHG supplementary figure dir |
| `figures/main/`, `figures/supplementary/` | New manuscript-numbered figure dirs |

---

## 2. Revised Figure 2 — k = 8 NMF component composition

### Passed

- File present: `figures/main/figure_2_nmf_k8_population_composition.{pdf,png}` ✓
- Source script self-describes as the k=8 replacement for the legacy k=2..5 branch figure ✓
- Population set: **EAS + EUR only**, 5 + 5 = 10 sub-populations (matches Lange 2013 Fig. 1A scope and the surrounding caption text) ✓
- Panel structure: (a) mean NMF loading stacked, (b) row-normalized composition stacked ✓
- Component labels in legend (script lines 47–56):
  - c4 = "Axis I  MICA" ✓
  - c5 = "Axis II  HLA-B / HLA-C" ✓
  - c6 = "Axis I  MICA" ✓
- Figure title (script line 204–207): "Figure 2.  k = 8 NMF component composition among rs2596542-T carrier haplotypes (EAS + EUR)" ✓
- Companion 26-pop figure exists for the supplement: `figures/supplementary/figure_2_nmf_k8_26pop.{pdf,png}` ✓
- Companion table exists: `tables/figure_2_component_population_summary.csv` (long format, 8 × 26 = 208 rows) ✓
- Figure-2 revision rationale documented at `reports/figure_2_revision_audit.md` ✓
- The text around the proposed Figure 2 caption (draft at `reports/figure_2_revision_audit.md` lines 79–81) says k = 8 consistently; no internal contradiction.

### Issues found

**ISSUE F2-1 — `docs/FIGURE_NUMBERING.md` row for Fig 2 still describes the OLD k=2..5 branch composition.**
Location: `docs/FIGURE_NUMBERING.md:15`.
Current row:

> `**Fig 2** | results/figure_2_branch_composition.{pdf,png} | cluster_anchor_haplotypes.py … | Anchor-carrier haplotype branch composition … k=2..5 hierarchical clustering …`

Recommended replacement (single row):

```markdown
| **Fig 2** | `figures/main/figure_2_nmf_k8_population_composition.{pdf,png}` | `figure_2_nmf_k8_population_composition.py` (consumes `results/nmf_W_26_k8.parquet` + `results/_supp_hap_meta.parquet`) | k = 8 NMF component composition among the 2,111 rs2596542-T carrier haplotypes, shown across the 5 EAS + 5 EUR sub-populations of Lange 2013 Fig. 1A. Panel (a) per-population mean NMF loading; panel (b) row-normalized component composition. **c4 + c6 = Axis I (MICA, EAS-enriched)**; **c5 = Axis II (HLA-B↑ / HLA-C↓, EUR-enriched)**. A 26-population companion is provided as Supplementary Figure 2; the legacy k = 2..5 hierarchical-clustering branch figure is archived as `figures/supplementary/figure_S_coarse_carrier_space_k2_k5.{pdf,png}` and is no longer cited in the main text. |
```

**ISSUE F2-2 — `README.md` (`the working repository/README.md`) still lists `figure_2_branch_composition.{pdf,png}` (line 167) and describes `cluster_anchor_haplotypes.py` as the figure-2 script (line 88).**
Recommended action: update the README repo-tree and per-script summary to point to `figure_2_nmf_k8_population_composition.py` and `figures/main/figure_2_nmf_k8_population_composition.{pdf,png}`. The branch-clustering script remains in the repo as a precursor diagnostic but should no longer be described as the Figure-2 source.

**ISSUE F2-3 — `docs/INTERNAL_NOTE.md` references k=2..5 branch composition in §"Method", §"Result", and the file-table at the bottom.**
This is an internal exploration note (marked "paper 1 shelved; this is intellectual closure"); the k=2..5 references are historically accurate within the note. **No action required** unless the note is published; in that case, add a header line explaining that branch composition is a precursor diagnostic, superseded in the manuscript by the k = 8 NMF view.

**ISSUE F2-4 — old Fig 2 PDF still present in `supplementary/figures/figure_2_branch_composition.{pdf,png}`.**
This file is the previous Figure 2 (k=2..5 hierarchical clustering) and is now a duplicate of the archived `figures/supplementary/figure_S_coarse_carrier_space_k2_k5.{pdf,png}` under a different filename. Recommended action: delete `supplementary/figures/figure_2_branch_composition.{pdf,png}` (or move to a `_archive/` subdirectory) to prevent ambiguity. The archived copy at `figures/supplementary/figure_S_coarse_carrier_space_k2_k5.{pdf,png}` is the canonical legacy version per `reports/figure_2_revision_audit.md`.

### No old "k = 2–5" wording in the new Figure 2 caption draft.

The draft caption in `reports/figure_2_revision_audit.md` says k = 8 throughout and only mentions k = 2..5 once, explicitly to explain that the lower-rank view is archived.

---

## 3. Supplementary Figure S7 — COMT Val158Met flip-flop replication / NMF extension

### Passed

- File present: `figures/supplementary/figure_S4_comt_flipflop.{pdf,png}` ✓
- Source script: `figure_S4_comt_flipflop.py`, reads from `data/comt_flipflop_nmf/tables/` ✓
- Three panels match the user-suggested final caption:
  - (A) Population-conditional signed LD (signed r) and asymmetry component C for rs4680 × rs2097603 across 1000G populations; blue = positive, brown = negative ✓
  - (B) Row-normalized k = 4 NMF component mixture conditional on rs4680 allele (V = Val/REF carrier, M = Met/ALT carrier) ✓
  - (C) Population-level COMT NMF mixture across 5,008 haplotypes ✓
- Suptitle in script (lines 234–237): "Supplementary Figure S7.  Reference-panel replication and NMF extension of the COMT Val158Met (rs4680) flip-flop example" ✓
- COMT is **not cited** in any of `docs/`, `supplementary/`, `reports/`, or `README.md` as a main-Results element (grep confirmed zero hits for "COMT" / "rs4680" / "Val158Met" / "Lin et al" outside the S7 figure source itself). This is consistent with the user's directive that COMT is proof-of-concept only ✓
- Supplemental Note S1 (if any) was not found as a separate file; the `supplementary/SUPPLEMENTARY_INDEX.md` and `docs/manuscript_revisions_v2.md` do not contain any "old S7" caption text that could leak into the new framing ✓

### Issues found

**ISSUE S7-1 — `supplementary/SUPPLEMENTARY_INDEX.md` lists only six supplementary figures (S1–S6); no S7 entry.**
Location: `supplementary/SUPPLEMENTARY_INDEX.md:9` ("Supplementary Figures (6)") and the file-system layout at lines 117–123.
Recommended action: add a new figure entry after the S6 block and update the count.

Suggested S7 entry (use the user's final caption verbatim):

```markdown
**Figure S7 — Reference-panel replication and NMF extension of the COMT Val158Met flip-flop example.**
Source: `figure_S4_comt_flipflop.py` (consumes `data/comt_flipflop_nmf/tables/comt_pair_metrics_by_pop.tsv`, `comt_component_mixture_by_rs4680_status.tsv`, `comt_component_mixture_by_pop.tsv`).
(A) Population-conditional signed LD and asymmetry component C for rs4680 (Val158Met) and rs2097603 across 1000 Genomes populations. Blue and brown indicate positive and negative signed values, respectively. (B) Row-normalized k = 4 NMF component mixture conditional on rs4680 allele status, with V indicating Val/REF carriers and M indicating Met/ALT carriers. (C) Population-level COMT NMF component mixture across all 5,008 haplotypes. Compared with the MHC rs2596542 locus, the COMT decomposition is simpler and more population-like, supporting broader applicability of carrier-set decomposition while illustrating locus-dependent differences in functional interpretability.
```

And in the file-system layout block (`supplementary/SUPPLEMENTARY_INDEX.md` lines 117–123):

```
supplementary/
├── SUPPLEMENTARY_INDEX.md
├── figures/
│   ├── figure_S1_lda_vs_nmf.{pdf,png}
│   ├── figure_S2_cophenetic_k.{pdf,png}
│   ├── figure_S3_seed_stability.{pdf,png}
│   ├── figure_S4_liri_weighted.{pdf,png}
│   ├── figure_S5_liri_hlaa_diagnostic.{pdf,png}
│   ├── figure_S6_liri_sensitivity.{pdf,png}
│   └── figure_S4_comt_flipflop.{pdf,png}
└── tables/
    └── ...
```

Also change the section header from `## Supplementary Figures (6)` to `## Supplementary Figures (7)`.

**ISSUE S7-2 — split supplementary figure directories.**
The pre-AJHG supplement uses `supplementary/figures/` (S1–S6 live there), but the S7 file is at `figures/supplementary/figure_S4_comt_flipflop.{pdf,png}` *and* a copy is now also under `supplementary/figures/figure_S4_comt_flipflop.{pdf,png}`. Two parallel directories with overlapping content is fragile.
Recommended action: pick one canonical location for the AJHG submission package — recommend `supplementary/figures/` (matches S1–S6 and `SUPPLEMENTARY_INDEX.md` layout block) — and have the source script write to that directory. Either remove `figures/supplementary/figure_S7_*` after the script is updated, or keep `figures/supplementary/` only for the new manuscript-numbered figures (Fig 2 k=8 supp, archived legacy). The two-directory split is the single biggest source of confusion in this revision.

**ISSUE S7-3 — S7 caption text not yet written into any caption file.**
The figure script embeds a suptitle (lines 234–237) and per-panel titles, but the figure caption itself does not appear in any of `docs/manuscript_revisions_v2.md`, `supplementary/SUPPLEMENTARY_INDEX.md`, or `reports/`. Apply the suggested caption (see S7-1) so the figure travels with its caption.

---

## 4. Supplementary Table S2 — top-signature SNVs

### Passed

- File present: `supplementary/tables/table_S2_component_top_snvs.csv` ✓
- Row count: 2840 = 8 components × 355 SNVs each ✓
- All 8 components are populated (c0..c7, 355 rows each) ✓
- 13 columns: `component, rank_within_component, position_grch37, H_loading, gene_region, C_flip, class_flip, eQTL_MICA, eQTL_MICB, eQTL_HLA-B, eQTL_HLA-C, eQTL_HCG27, eQTL_HCP5` ✓
- The anchor SNV rs2596542 (chr6:31,366,595, GRCh37) is **absent** from the table (0 rows match), confirming anchor exclusion ✓
- `gene_region` takes a clean closed vocabulary: `MICA / MICB / HLA-B / HLA-C / HCG27 / HCP5 / 3'_of_MICB / 5'_of_HCG27 / intergenic_class_I` ✓
- Top-5% / 95th-percentile definition is documented in `docs/manuscript_revisions_v2.md` §B3 (Methods text replacement) — re-quoted in the suggested table legend below ✓
- Anchor removal is documented in `docs/manuscript_revisions_v2.md` §B5 (Methods text) ✓

### Issues found

**ISSUE S2-1 — `H_loading` is not rounded.**
Sample rows show full IEEE-754 precision (`0.9388098120689392`, `0.9348764419555664`). The user's recommended legend says "H-matrix loading values are rounded reasonably". Two options:
- (a) Round `H_loading` to 4 decimal places in `supplementary/tables/table_S2_component_top_snvs.csv` (changes 2,840 cells; cheap one-liner in pandas or awk).
- (b) Drop "rounded reasonably" from the legend and report values at full precision so they're byte-identical with the parquet H matrix.

Recommended option: (a), to match the legend.

Patch (one-liner, run from `mica_flipflop/`):

```bash
python3 -c "
import pandas as pd
p = 'supplementary/tables/table_S2_component_top_snvs.csv'
df = pd.read_csv(p)
df['H_loading'] = df['H_loading'].round(4)
df.to_csv(p, index=False)
print('rounded H_loading to 4 decimals; rows=', len(df))
"
```

**ISSUE S2-2 — `SUPPLEMENTARY_INDEX.md` Table S2 description states "ordered by H_loading" but the file shows ordering by `rank_within_component` (rank ascending within component).**
This is consistent in practice (rank is monotonic with descending H_loading), but the wording "ordered by H_loading" could be misread as descending-by-loading. Recommend the legend explicitly say "ordered by descending H_loading within each component" (rank 1 = highest loading).

**ISSUE S2-3 — Table S2 legend should be appended to `SUPPLEMENTARY_INDEX.md`'s Table S2 paragraph.**

Recommended Table S2 legend (final form, drop-in to `SUPPLEMENTARY_INDEX.md:40-41`):

> **Table S2 — Top-signature SNVs defining NMF components in rs2596542-T carrier haplotypes.**
> Top-signature SNVs were defined as variants with H-matrix loadings in the top 5% of each component's H[c, ·] distribution (i.e., at or above the 95th percentile), yielding 355 SNVs per component. The table lists the component-defining SNV sets used for LD-reversal enrichment, *cis*-eQTL annotation, MAF×distance matched-permutation analysis (Supplementary Table S3), and tag-SNP component mapping (Supplementary Table S8). The anchor SNV rs2596542 was removed from the NMF input matrix after carrier selection because all retained haplotypes carry the rs2596542-T allele at this site, leaving it with zero variance in the analytic matrix. Columns: `component` (c0..c7); `rank_within_component` (1 = highest H_loading); `position_grch37` (chr6 1000G phase 3 coordinate); `H_loading` (rounded to 4 decimal places); `gene_region` (MICA / MICB / HLA-B / HLA-C / HCG27 / HCP5 / 3'_of_MICB / 5'_of_HCG27 / intergenic_class_I); `C_flip`, `class_flip` (boolean LD-reversal flags); `eQTL_<gene>` (boolean indicating whether the SNV is a GTEx Whole_Blood *cis*-eQTL for that gene). Rows: 2,840 (8 components × 355).

**ISSUE S2-4 — Whether Table S2 should be pasted into the Word manuscript.**
2,840 rows × 13 columns ≈ 37,000 cells. **Do not paste into Word.** Reference the CSV in the supplementary materials; provide the CSV as an upload to the journal's supplementary file system. Add a one-line statement in the manuscript Methods: "The full top-signature SNV table (8 components × 355 SNVs) is provided as Supplementary Table S2 (`table_S2_component_top_snvs.csv`)."

---

## 5. Numbering and cross-reference audit

### Supplementary figures (target: S1–S7)

| # | Expected | File present | Status |
|---|---|---|---|
| S1 | NMF vs LDA | `supplementary/figures/figure_S1_lda_vs_nmf.{pdf,png}` | ✓ |
| S2 | NMF rank selection | `supplementary/figures/figure_S2_cophenetic_k.{pdf,png}` | ✓ |
| S3 | NMF seed stability | `supplementary/figures/figure_S3_seed_stability.{pdf,png}` | ✓ |
| S4 | LIRI-JP c5 weighting sensitivity | `supplementary/figures/figure_S4_liri_weighted.{pdf,png}` | ✓ |
| S5 | HLA-A adjustment diagnostic | `supplementary/figures/figure_S5_liri_hlaa_diagnostic.{pdf,png}` | ✓ |
| S6 | LIRI-JP regression sensitivity + multiple-testing | `supplementary/figures/figure_S6_liri_sensitivity.{pdf,png}` | ✓ |
| S7 | COMT replication / NMF extension | `figures/supplementary/figure_S4_comt_flipflop.{pdf,png}` + duplicate copy at `supplementary/figures/figure_S4_comt_flipflop.{pdf,png}` | ✓ (location split — see S7-2) |

### Duplicate-S6 question

Two files with overlapping S6 semantics exist in `supplementary/figures/`:

- `figure_S6_liri_sensitivity.{pdf,png}` ← canonical S6 (matches `SUPPLEMENTARY_INDEX.md`)
- `figure_supp_6_liri_sensitivity.{pdf,png}` ← legacy filename, byte-identical or near-identical content

**ISSUE N-1 — delete the legacy `figure_supp_6_liri_sensitivity.{pdf,png}` copies** to remove ambiguity; the canonical S6 filename is `figure_S6_liri_sensitivity.{pdf,png}`.

### Supplementary tables (target: S1–S10)

| # | Expected | File(s) | Status |
|---|---|---|---|
| S1 | LIRI-JP regression | `table_S1_regression_models.csv` + `table_S1d_VIF.csv` | ✓ |
| S2 | Top-signature SNVs | `table_S2_component_top_snvs.csv` | ✓ |
| S3 | Matched permutation | `table_S3_matched_permutation.csv` | ✓ |
| S4 | GTEx tissue eQTL NES | `table_S4_gtex_tissue_eqtl_NES.csv` | ✓ |
| S5 | LDA vs NMF | `table_S5_lda_vs_nmf.csv` | ✓ |
| S6 | NMF rank selection | `table_S6_nmf_rank_cophenetic.csv` | ✓ |
| S7 | NMF seed stability | `table_S7_nmf_seed_stability.csv` + `table_S7b_nmf_seed_stability_per_seed.csv` | ✓ |
| S8 | Tag SNP component mapping | `table_S8_tag_snp_component_mapping.csv` + `_wide.csv` | ✓ |
| S9 | LIRI sensitivity + multiple testing | `table_S9_liri_sensitivity.csv` + `table_S9b_multiple_testing.csv` | ✓ |
| S10 | Axis sign coherence | `table_S10_axis_sign_coherence.csv` | ✓ |

No missing supplementary table number. **The supplement now contains seven figures (S1–S7) and ten tables (S1–S10), not six and nine as `SUPPLEMENTARY_INDEX.md` headers still suggest.**

### Figure-3 cross-reference

`FIGURE_NUMBERING.md:19` ("Fig 3 update (panel c)") correctly cross-references "Supplementary Table S10 (`table_S10_axis_sign_coherence.csv`)" — consistent with the file present. No fix needed.

### Main-Results cross-references

`docs/manuscript_revisions_v2.md` cites **Tables S1, S3, S4, S5, S6, S7, S8, S9, S9b, S10** and **Figures S1, S2, S3, S4, S5, S6**. Notably:

- It does **not** yet cite Supplementary Figure S7 anywhere — consistent with the user's directive that COMT is proof-of-concept and not cited in main Results.
- It does cite Table S2 implicitly via Methods §B3 (top-5% definition) — sufficient. An explicit "Supplementary Table S2" cross-reference should be added to the NMF Methods sub-section as: "...yielding 355 SNVs per component (Supplementary Table S2)."

### Results-text references claim from the task brief

> Results text references Supplemental Figs. S1–S6 and Tables S1 / S7 / S9 / S10 consistently

This is true in `docs/manuscript_revisions_v2.md`. The new Figure S7 (COMT) is **not** cited from Results, by design.

---

## 6. Summary table

| Item | Pass | Issue | Fix priority |
|---|---|---|---|
| Fig 2 file in `figures/main/` | ✓ | — | — |
| Fig 2 EAS+EUR + 26-pop supp | ✓ | — | — |
| Fig 2 c4/c6/c5 labels | ✓ | — | — |
| Fig 2 k=8 wording | ✓ | — | — |
| `FIGURE_NUMBERING.md` Fig 2 row | ✗ | F2-1 | high |
| `README.md` Fig 2 references | ✗ | F2-2 | high |
| Legacy `figure_2_branch_composition.pdf` still in `supplementary/figures/` | ✗ | F2-4 | medium (delete or move) |
| Fig S7 file + caption text | ✓ | — | — |
| Fig S7 in `SUPPLEMENTARY_INDEX.md` | ✗ | S7-1 | high |
| Fig S7 single canonical location | ✗ | S7-2 | medium |
| Fig S7 caption written in caption file | ✗ | S7-3 | high (fixed by applying S7-1) |
| COMT not cited as main-Results locus | ✓ | — | — |
| Table S2 file structure | ✓ | — | — |
| Table S2 anchor removal | ✓ | — | — |
| Table S2 H_loading rounding | ✗ | S2-1 | medium (round to 4 d.p.) |
| Table S2 legend wording | ✗ | S2-2, S2-3 | high |
| Figs S1–S6, S7 present | ✓ | — | — |
| Tables S1–S10 present | ✓ | — | — |
| Duplicate S6 (`figure_supp_6_*`) | ✗ | N-1 | medium (delete) |
| `SUPPLEMENTARY_INDEX.md` headers say "6 figures / 9 tables" | ✗ | (within S7-1) | high |

---

## 7. Unresolved items (require user input)

- **U-1 — Where is the manuscript DOCX?** The git-tracked content in this repo (`docs/`, `supplementary/`, `figures/`, source scripts) does not contain a manuscript Word file or a master `manuscript.md` / `main.tex`. `docs/manuscript_revisions_v2.md` is a *patch set* (search-and-replace recipe), not the manuscript itself. The full audit therefore cannot verify the exact wording in the manuscript; this add-on audit verifies only the repo-side support materials. **Confirm the location of the manuscript file (Word, LaTeX, or markdown) and whether it should be added to this repo or remain external.**
- **U-2 — Should the archived legacy Fig 2 (`figure_S_coarse_carrier_space_k2_k5`) be promoted to a numbered slot (e.g., Supplementary Figure S8)?** Currently it has an unnumbered descriptive filename in `figures/supplementary/`. Decide before submission.
- **U-3 — Should `supplementary/figures/` and `figures/supplementary/` be consolidated?** See S7-2. Strong recommendation: yes, into `supplementary/figures/` to match `SUPPLEMENTARY_INDEX.md`.
- **U-4 — Is `tables/figure_2_component_population_summary.csv` (208-row companion to revised Fig 2) intended as a supplementary table?** If yes, suggest assigning a number (e.g., Supplementary Table S11) and adding a legend; if no, leave as a working artifact alongside the figure source.

---

## 8. Recommended fix order (smallest blast radius first)

1. Update `supplementary/SUPPLEMENTARY_INDEX.md` headers and add S7 entry. (Issue S7-1.)
2. Round `H_loading` in `table_S2_component_top_snvs.csv` to 4 d.p., and update the Table S2 legend with the wording in §4 above. (Issues S2-1, S2-2, S2-3.)
3. Delete the duplicate `supplementary/figures/figure_supp_6_liri_sensitivity.{pdf,png}`. (Issue N-1.)
4. Decide canonical location of S7 (`supplementary/figures/`), update `figure_S4_comt_flipflop.py` to write there, remove the duplicate copy in `figures/supplementary/`. (Issue S7-2.)
5. Update `docs/FIGURE_NUMBERING.md` Fig 2 row and the README repo-tree to point to `figures/main/figure_2_nmf_k8_population_composition.{pdf,png}`. (Issues F2-1, F2-2.)
6. Delete (or move to `_archive/`) the old `supplementary/figures/figure_2_branch_composition.{pdf,png}`. (Issue F2-4.)
7. (After U-2 decision) Renumber `figure_S_coarse_carrier_space_k2_k5` if it is to be cited as a numbered supplementary figure.

---

End of add-on audit.
