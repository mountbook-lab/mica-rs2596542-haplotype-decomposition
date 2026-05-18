# GitHub update summary — final supplementary figure/table revisions

Date: 2026-05-15 (audit) → 2026-05-15 (fixes applied)
Repo (working): `~/paper1-proxy-topology/mica_flipflop/`
Repo (GitHub mirror): `~/github_repos/mica-rs2596542-haplotype-decomposition/` → `mountbook-lab/mica-rs2596542-haplotype-decomposition`
Last GitHub commit: `2c82980 Add Zenodo DOI to README` (2026-05-08)

This document lists what should be brought into the GitHub mirror to
publish the post-2026-05-08 manuscript-revision support materials. It
is not a `git diff` because `~/paper1-proxy-topology/` is currently
**uninitialised for commits** (no commits yet on `master`); changes
must be brought to the GitHub clone manually or by `rsync`/`cp`, then
committed there.

## Applied fixes (2026-05-15)

The following audit issues from `reports/addon_audit_final_figures_tables.md`
have been resolved in the working repo. No git operations performed on
either repo; cross-copy + commit remains a manual step.

| Audit ID | Fix | File(s) |
|---|---|---|
| **F2-1** | FIGURE_NUMBERING.md Fig 2 row replaced with k = 8 NMF wording | `docs/FIGURE_NUMBERING.md` |
| **F2-2** | README.md repo-tree + script table updated to point to `figures/main/figure_2_nmf_k8_population_composition.{pdf,png}` | `README.md` |
| **F2-4** | Old `figure_2_branch_composition.{pdf,png}` moved to `supplementary/figures/_archive/` | filesystem |
| **S7-1** | SUPPLEMENTARY_INDEX.md header counts updated to "7 figures, 10 tables"; new S7 entry added with the user-supplied caption | `supplementary/SUPPLEMENTARY_INDEX.md` |
| **S7-2** | S7 PDF/PNG consolidated into `supplementary/figures/` (canonical); duplicate copies removed from `figures/supplementary/`; `figure_S7_comt_flipflop.py` OUT_DIR updated | filesystem + `figure_S7_comt_flipflop.py` |
| **S7-3** | S7 caption written into SUPPLEMENTARY_INDEX.md | `supplementary/SUPPLEMENTARY_INDEX.md` |
| **N-1** | Legacy duplicate `figure_supp_6_liri_sensitivity.{pdf,png}` deleted from `supplementary/figures/` | filesystem |
| **S2-1** | Table S2 `H_loading` rounded to 4 decimal places (preserves all 2,840 rows and 13 columns) | `supplementary/tables/table_S2_component_top_snvs.csv` |
| **S2-2, S2-3** | Table S2 legend rewritten with rounding note, anchor-removal rationale, and full column dictionary; `make_supplementary_package.py` updated to emit the new legend and to round on write | `supplementary/SUPPLEMENTARY_INDEX.md` + `make_supplementary_package.py` |
| **(new)** | `make_supplementary_package.py` "Supplementary Figures (6)" / "Supplementary Tables (9)" string literals updated to "(7)" / "(10)" so re-runs do not revert the INDEX | `make_supplementary_package.py` |
| **(new)** | Legacy `figure_S_coarse_carrier_space_k2_k5.{pdf,png}` also moved to `supplementary/figures/_archive/` (user decision: no number assigned) | filesystem |
| **Inventory** | `manuscript/figure_table_inventory_final.md` "Outstanding decisions" section updated to record U-2/U-3/U-4 resolutions | `manuscript/figure_table_inventory_final.md` |

### Items intentionally left for manual follow-up

- `make_supplementary_package.py` still emits only S1–S6 figures and S1–S9 tables in its INDEX template body (only the count strings were patched). The current `supplementary/SUPPLEMENTARY_INDEX.md` has been manually extended with the S7 + S10 paragraphs. **If the generator is re-run, the S7 figure and S10 table entries will need to be re-added manually**, or the generator body extended with S7 / S10 templates.
- `supp_6_liri_sensitivity.py:34, 341` still writes its render to `figure_supp_6_liri_sensitivity.{pdf,png}`; `make_supplementary_package.py:478` then renames to the canonical `figure_S6_liri_sensitivity.{pdf,png}`. Acceptable two-step pipeline — left as-is.
- `cluster_anchor_haplotypes.py` is preserved as a precursor diagnostic. Its README row now flags it as no longer the Figure 2 source.
- `docs/INTERNAL_NOTE.md` references to k = 2..5 branch composition remain — historically accurate within the internal note; not for publication.

---

## 1. New files to bring to the GitHub clone

### 1.1 New source scripts

| File | Purpose | Priority |
|---|---|---|
| `figure_2_nmf_k8_population_composition.py` | Source for the revised Figure 2 | **must** |
| `figure_S7_comt_flipflop.py` | Source for Supplementary Figure S7 (COMT) | **must** |
| `figure_style.py` | Shared AJHG style helpers used by the above | **must** |
| `supp_common.py` | Shared loader used by `supp_1..6` and Fig 2 | **must** |
| `supp_1_lda_comparison.py` | Supp Fig S1 / Table S5 source | **must** |
| `supp_2_cophenetic_k.py` | Supp Fig S2 / Table S6 source | **must** |
| `supp_3_seed_stability.py` | Supp Fig S3 / Table S7 source | **must** |
| `supp_4_liri_weighted.py` | Supp Fig S4 / Table S1B source | **must** |
| `supp_5_liri_hlaa_conditioning.py` | Supp Fig S5 source | **must** |
| `supp_6_liri_sensitivity.py` | Supp Fig S6 / Table S9 source | **must** |
| `make_supplementary_package.py` | Builds Tables S1, S2, S8 from raw results | **must** |
| `axis_sign_coherence_matrix.py` | Fig 3 panel (c) / Table S10 source | **must** |
| `axis_I_cross_tissue.py` | Cross-tissue NES check feeding Table S10 | **must** |
| `tag_snp_component_map.py` | Table S8 source (already present in github_repos under same name — verify identical) | should compare |

### 1.2 New manuscript-numbered figure files

| Path | Description |
|---|---|
| `figures/main/figure_2_nmf_k8_population_composition.pdf` | Revised main Figure 2 |
| `figures/main/figure_2_nmf_k8_population_composition.png` | (raster) |
| `figures/supplementary/figure_2_nmf_k8_26pop.pdf` | 26-pop companion |
| `figures/supplementary/figure_2_nmf_k8_26pop.png` | (raster) |
| `figures/supplementary/figure_S_coarse_carrier_space_k2_k5.pdf` | Archived legacy Fig 2 |
| `figures/supplementary/figure_S_coarse_carrier_space_k2_k5.png` | (raster) |
| `supplementary/figures/figure_S1_lda_vs_nmf.{pdf,png}` | Supp Fig S1 |
| `supplementary/figures/figure_S2_cophenetic_k.{pdf,png}` | Supp Fig S2 |
| `supplementary/figures/figure_S3_seed_stability.{pdf,png}` | Supp Fig S3 |
| `supplementary/figures/figure_S4_liri_weighted.{pdf,png}` | Supp Fig S4 |
| `supplementary/figures/figure_S5_liri_hlaa_diagnostic.{pdf,png}` | Supp Fig S5 |
| `supplementary/figures/figure_S6_liri_sensitivity.{pdf,png}` | Supp Fig S6 |
| `supplementary/figures/figure_S7_comt_flipflop.{pdf,png}` *or* `figures/supplementary/figure_S7_comt_flipflop.{pdf,png}` | Supp Fig S7 (canonical location to be decided — see audit issue S7-2) |

### 1.3 New supplementary tables

| Path |
|---|
| `supplementary/tables/table_S1_regression_models.csv` |
| `supplementary/tables/table_S1d_VIF.csv` |
| `supplementary/tables/table_S2_component_top_snvs.csv` |
| `supplementary/tables/table_S3_matched_permutation.csv` |
| `supplementary/tables/table_S4_gtex_tissue_eqtl_NES.csv` |
| `supplementary/tables/table_S5_lda_vs_nmf.csv` |
| `supplementary/tables/table_S6_nmf_rank_cophenetic.csv` |
| `supplementary/tables/table_S7_nmf_seed_stability.csv` |
| `supplementary/tables/table_S7b_nmf_seed_stability_per_seed.csv` |
| `supplementary/tables/table_S8_tag_snp_component_mapping.csv` |
| `supplementary/tables/table_S8_tag_snp_component_mapping_wide.csv` |
| `supplementary/tables/table_S9_liri_sensitivity.csv` |
| `supplementary/tables/table_S9b_multiple_testing.csv` |
| `supplementary/tables/table_S10_axis_sign_coherence.csv` |
| `tables/figure_2_component_population_summary.csv` |

### 1.4 New documentation

| Path | Purpose |
|---|---|
| `supplementary/SUPPLEMENTARY_INDEX.md` | Reader-facing supplementary index (apply audit issue S7-1 first) |
| `docs/FIGURE_NUMBERING.md` | Locked figure-numbering map (apply audit issue F2-1 first) |
| `docs/manuscript_revisions_v2.md` | Drop-in manuscript revision recipes |
| `docs/rs2596542_allele_strand_build_mapping.md` | Allele / strand / build documentation |
| `docs/rs2596542_allele_strand_build_mapping.csv` | Companion CSV |
| `reports/figure_2_revision_audit.md` | Documents the k=2..5 → k=8 figure switch |
| `reports/addon_audit_final_figures_tables.md` | This add-on audit |
| `reports/github_update_summary_final.md` | This file |

### 1.5 New figure/table inventory

| Path |
|---|
| `manuscript/figure_table_inventory_final.md` (new file; see this audit's companion deliverable) |

---

## 2. Files to **exclude** from GitHub

These are work-in-progress, paper-2-scoped, or contain large or
controlled-access data and **should not** be pushed to the public repo.

| Path | Reason |
|---|---|
| `_future_liri_paper/` | Paper-2 scope; deliberately quarantined per project notes ([[mica_flipflop_liri_supp]] memory) |
| `data/region_500kb.vcf.gz`, `data/region_500kb.vcf.gz.tbi`, `data/three_snps.vcf` | Derived from 1000 Genomes (publicly redownloadable); not the canonical distribution channel |
| `results/` (112 MB) | Working caches (parquet / CSV intermediates). Keep `results/*_summary.csv` if cited from text; otherwise skip the directory in favour of the published `supplementary/tables/` |
| `__pycache__/`, `*.pyc` | Build artifacts |
| `paper2_supp_section3_baseline.py`, `paper2_supp_section4a_component_by_pop.py`, `paper2_supp_section4b_k6_internal_structure.py` | Paper-2 scope |
| `c3_tissue_micb_check.py`, `liri_tumor_vs_normal_c4.py`, `gse223201_lenvatinib_axis.py`, `ins_anchor_*.py`, `explore_rs11509487_*.py`, `axis_I_cross_tissue.py` (if not feeding a published figure) | Exploratory / paper-2 reuse; review individually |
| `docs/INTERNAL_NOTE.md` | Internal exploration; not for publication (marked "paper 1 shelved") |
| `supplementary/figures/figure_supp_6_liri_sensitivity.{pdf,png}` | Legacy duplicate of canonical S6 (audit issue N-1) |
| `supplementary/figures/figure_2_branch_composition.{pdf,png}` | Legacy main Fig 2 in old location (audit issue F2-4) |
| `supplementary/figures/figure_1_lange_replication.{pdf,png}` … `figure_5_liri_panel.{pdf,png}` | These are main figures duplicated into the supp dir; prefer `results/figure_*.{pdf,png}` (the source-of-truth main figure location per `FIGURE_NUMBERING.md`) |
| One of the two S7 copies | See audit issue S7-2; keep one canonical copy only |

---

## 3. Modifications status — all applied 2026-05-15

All pre-existing files flagged in the original audit have been brought
in line with the revised Fig 2 / new S7 / new S10. See "Applied fixes
(2026-05-15)" table at the top of this document. Deletions and file
moves are summarised in §2 and were performed on the working repo.

---

## 4. Suggested commit sequence on the GitHub mirror

To keep history readable, split the work into three commits rather
than one mega-commit:

### Commit 1 — manuscript revision deliverables

```
Finalize supplementary figure/table updates and add targeted audit

- Replace main Figure 2 with k = 8 NMF component composition
  (figures/main/figure_2_nmf_k8_population_composition.{pdf,png}).
  Source: figure_2_nmf_k8_population_composition.py. Companion 26-pop
  version under figures/supplementary/figure_2_nmf_k8_26pop.{pdf,png}.
  The previous k = 2..5 hierarchical-clustering figure is archived as
  figures/supplementary/figure_S_coarse_carrier_space_k2_k5.{pdf,png}.

- Add Supplementary Figure S7 (COMT Val158Met flip-flop reference-
  panel replication and NMF extension; figure_S7_comt_flipflop.py).
  COMT is an external proof-of-concept; not cited from main Results.

- Add supplementary support scripts (supp_1..6, supp_common,
  axis_sign_coherence_matrix, axis_I_cross_tissue,
  make_supplementary_package) and the resulting tables/figures
  (Supplementary Tables S1–S10, Supplementary Figures S1–S6).

- Update supplementary/SUPPLEMENTARY_INDEX.md to "7 figures, 10
  tables" and add the S7 entry.

- Update docs/FIGURE_NUMBERING.md Fig 2 row and README.md repo-tree
  to point to the revised Figure 2.

- Round table_S2_component_top_snvs.csv H_loading to 4 d.p. and
  expand the Table S2 legend.

- Add reports/figure_2_revision_audit.md,
  reports/addon_audit_final_figures_tables.md,
  reports/github_update_summary_final.md, and
  manuscript/figure_table_inventory_final.md.
```

### Commit 2 — cleanup of legacy/duplicate files

```
Remove legacy and duplicate supplementary figure files

- Delete supplementary/figures/figure_supp_6_liri_sensitivity.{pdf,png}
  (legacy duplicate of figure_S6_liri_sensitivity.{pdf,png}).
- Delete supplementary/figures/figure_2_branch_composition.{pdf,png}
  (legacy main Figure 2; archived copy retained at
  figures/supplementary/figure_S_coarse_carrier_space_k2_k5.{pdf,png}).
- Consolidate Supplementary Figure S7 to a single canonical location
  in supplementary/figures/figure_S7_comt_flipflop.{pdf,png}.
```

### Commit 3 — documentation alignment (optional)

```
Sync README and FIGURE_NUMBERING references to revised Figure 2
```

If the user prefers a single commit, use commit-1's message verbatim as
the headline; the suggested commit message from the task brief
("Finalize supplementary figure/table updates and add targeted audit")
matches this exactly.

---

## 5. Pre-push checklist

- [ ] Apply the five `modify before commit` items in §3.
- [ ] Run audit-recommended one-liner to round Table S2 `H_loading` (see audit §4, issue S2-1).
- [ ] Verify `figure_2_nmf_k8_population_composition.py` runs end-to-end from a clean checkout (its inputs `results/nmf_W_26_k8.parquet` and `results/_supp_hap_meta.parquet` are not committed — either include the parquets as Git LFS / release artifacts, or document the upstream `decompose_branches_26.py` step in the README).
- [ ] Same for `figure_S7_comt_flipflop.py` — its inputs live in `~/analyses/comt_flipflop_nmf/tables/`. Either commit those TSVs alongside the script, or rewrite the script to bundle them under `data/comt_flipflop_nmf/`.
- [ ] Check `figures/supplementary/figure_S7_comt_flipflop.{pdf,png}` is the freshly rendered version (matches script in this repo, not a stale older render).
- [ ] Confirm `.gitignore` excludes `__pycache__/` and `results/` if those are not to be committed.
- [ ] Do **not** push without explicit user authorization. The task brief explicitly says: *"Do not push automatically."*

---

End of GitHub update summary.
