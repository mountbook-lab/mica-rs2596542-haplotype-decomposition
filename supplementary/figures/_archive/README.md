# Archived figures — not cited from the manuscript

Supersedes the earlier legacy k = 2..5 branch-composition Figure 2; the
current Figure 2 is the k = 8 NMF component-composition view, and the
legacy branch-composition render is retained only in the archive for
reproducibility.

## Contents

| File | Description |
|---|---|
| `figure_2_branch_composition.{pdf,png}` | Original main Figure 2 from earlier drafts. Anchor-carrier (rs2596542-T) hierarchical clustering on Hamming distance, k = 2..5 cut points, rendered by `cluster_anchor_haplotypes.py`. |
| `figure_S_coarse_carrier_space_k2_k5.{pdf,png}` | Renamed copy of the same legacy figure produced during the k = 2..5 → k = 8 transition. Kept alongside the original so that prior drafts and rebuttal materials referencing either filename remain resolvable. |

## Why these are kept

These renders are no longer cited in the manuscript, but they document
the earlier visualisation choice and let reviewers or follow-up authors
reproduce the precursor diagnostic without having to re-run
`cluster_anchor_haplotypes.py` from raw data. The current Figure 2
(k = 8 NMF component composition) is at
`figures/main/figure_2_nmf_k8_population_composition.{pdf,png}` and is
described in `docs/FIGURE_NUMBERING.md`.

See also `reports/figure_2_revision_audit.md` for the full record of
the figure switch.
