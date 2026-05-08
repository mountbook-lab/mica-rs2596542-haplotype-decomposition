# Internal exploration note — does Lange 2013's signed-r flip decompose into a topology-class flip?

Date: 2026-05-08
Status: internal exploration only (paper 1 shelved; this is intellectual closure)

## Question

Lange 2013 (J Hepatol) showed that the signed Pearson correlation r
between rs2596542 (MICA upstream) and two nearby SNPs flips sign
between Japanese and European populations:

| pair | r in JPT | r in CEU |
|---|---|---|
| rs2596542 × rs2244546 | +0.59 | −0.19 |
| rs2596542 × rs9275572 | +0.27 | −0.12 |

Does this signed-r flip correspond to a *carrier-set topology class*
flip (DISJOINT / STRICT / IDENTICAL / PARTIAL), or is it independent
information that signed r captures and class doesn't?

If carrier-set topology adds new information, the v11 / paper-1
framework is conceptually richer than signed correlation alone.
If not, the framework just re-encodes signed correlation in different
words.

## Method

For each pair, computed across 10 EAS + EUR 1000G phase 3
sub-populations (CHB, JPT, CHS, CDX, KHV, CEU, TSI, FIN, GBR, IBS)
plus EAS_pooled / EUR_pooled:

- n11, n10, n01, n00 (haplotype-level co-occurrence counts)
- p_a, p_b (allele frequencies)
- signed Pearson r (REF/ALT = dbSNP forward; ALT coded as 1)
- P(b|a) = n11/(n11+n10), P(a|b) = n11/(n11+n01)
- Δ = P(b|a) − P(a|b) (asymmetry index)
- topology class: DISJOINT / STRICT_A_IN_B / STRICT_B_IN_A / IDENTICAL / PARTIAL

Anchor = rs2596542. Partners as above.

## Result — mixed answer

### Pair 1: rs2596542 × rs2244546 — class change DOES capture the flip

| pop | r_signed | n11 | class |
|---|---|---|---|
| JPT | +0.489 | 29 | PARTIAL |
| CHB | +0.067 | 6 | PARTIAL |
| CHS | +0.090 | 8 | PARTIAL |
| KHV | +0.108 | 8 | PARTIAL |
| CDX | +0.064 | 8 | PARTIAL |
| **CEU** | **−0.187** | **0** | **DISJOINT** |
| TSI | −0.143 | 3 | PARTIAL |
| FIN | −0.202 | 4 | PARTIAL |
| GBR | −0.163 | 1 | PARTIAL |
| IBS | −0.113 | 2 | PARTIAL |

- CEU has n11 = 0 → **DISJOINT** (not just PARTIAL with negative r).
- Other EUR populations have n11 = 1–4 → PARTIAL but very near the
  DISJOINT boundary.
- The signed-r flip from positive (EAS) to negative (EUR) is
  *concomitant* with a transition into the DISJOINT regime.
- Lange's signed r: +0.59 JPT (1000G: +0.49), −0.19 CEU (1000G: −0.19).
  Magnitudes are slightly different (HapMap3 vs 1000G phase 3 cohort
  composition), signs match.

**Topology class adds information here**: the CEU n11=0 is a true
DISJOINT class assignment, not just "small positive correlation
flipped to small negative". The r alone says "flipped"; the class
says "flipped *into* the regime where the carrier sets no longer
co-occur at all".

### Pair 2: rs2596542 × rs9275572 — class change does NOT capture the flip

| pop | r_signed | n11 | class |
|---|---|---|---|
| CHB | −0.274 | 27 | PARTIAL |
| JPT | −0.483 | 23 | PARTIAL |
| CHS | −0.324 | 26 | PARTIAL |
| CDX | −0.384 | 25 | PARTIAL |
| KHV | −0.151 | 41 | PARTIAL |
| CEU | +0.260 | 47 | PARTIAL |
| TSI | +0.115 | 60 | PARTIAL |
| FIN | +0.278 | 54 | PARTIAL |
| GBR | +0.251 | 58 | PARTIAL |
| IBS | +0.126 | 69 | PARTIAL |

- All 10 populations sit firmly in PARTIAL (every cell has
  n11, n10, n01, n00 substantially > 0).
- Sign flip in r is real (EAS negative, EUR positive in dbSNP-forward
  ALT-as-1 coding) but **not visible at the class level**.
- The flip is captured by signed r, by Δ shape changes
  (P(b|a) and P(a|b) values shift), but not by a class transition.

(Note on Lange's signs vs ours: Lange used "minor allele" coding;
rs9275572 has G major in both EAS and EUR, so my dbSNP-ALT=G coding
is opposite of Lange's, and signs flip. Re-coding to Lange's
convention reproduces +0.48 JPT, −0.26 CEU — same flip pattern,
moderately different magnitude than Lange's HapMap3 values.)

## Interpretation

**Carrier-set topology class is *partially* richer than signed r,
but not always.**

- When the signed-r flip is accompanied by *vanishing n11* (one of the
  populations effectively crosses into the DISJOINT regime), the class
  change is the more rigorous statement and captures something signed
  r merely flags. Pair 1 is this case.
- When the signed-r flip happens *within* the PARTIAL regime
  (all populations retain substantial co-occurrence), the topology
  class is constant and provides no extra information; only signed r
  (or the Δ asymmetry) distinguishes the populations. Pair 2 is this
  case.

So the framework's added value is **regime-dependent**: it shines when
populations approach a boundary (DISJOINT, IDENTICAL, STRICT
inclusion), and reduces to signed r when they don't.

This matches what paper 1's anchor-fixed framework was actually
demonstrating — boundary-class membership is the key analytic claim,
not the broader topology label. In the *boundary regime* the
framework adds information; in the *partial-overlap regime* signed
correlation already captures it.

## What this means for the (now-shelved) project

- Lange 2013 is correct that signed correlation flip is the headline
  observation; it captures the population heterogeneity in 9 out of
  10 cases tested here.
- Paper 1's framework is a useful refinement *only at the boundary
  regime* (DISJOINT, STRICT inclusion) — which is exactly where the
  paper's case studies sit (MICB rs11509487, etc).
- For a hypothetical reframe of the project: the genuine novelty is
  *boundary-regime tagging* — when does the tag SNP land in a regime
  where the carrier-set topology forces a binary class assignment, and
  what does that imply for cross-population transferability of
  associations?
- For Lange's example (rs2244546 in CEU n11=0), this is exactly the
  rare-anchor-blind-spot of paper 1 §7.3 (CEU rs2596542-T frequency
  ≈ 0.30, but rs2244546 G frequency ≈ 0.08, n11 expected ≈ 24 under
  independence vs observed 0 — a true DISJOINT signal).

## Files

| File | Purpose |
|---|---|
| `data/three_snps.vcf` | rs2596542, rs2244546, rs9275572 raw 1000G genotypes (chr6:31,366,595, 31,435,833, 32,678,999 GRCh37) |
| `analyze_topology_classes.py` | computes per-pop n11/n10/n01, signed r, P(v\|tag), P(tag\|v), Δ, class |
| `render_figure.py` | 2-panel figure: signed r vs topology class |
| `results/per_pop_topology.csv` | structured output |
| `results/figure_signed_r_vs_topology.{pdf,png}` | rendered figure |
| `INTERNAL_NOTE.md` | this file |

## Extension: region-wide M+C scan (added 2026-05-08)

The 2-SNP picture was followed up with a full anchor-fixed scan of
**every biallelic SNV in chr6:31,116,595–31,616,595 (±250 kb around
rs2596542)** — n = 20,019 partners, of which 10,000 had data in both
EAS and EUR sub-populations.

### Δ = M + C decomposition

- M = p_A − p_B (marginal frequency-difference component)
- C = Δ − M = D · (1/p_B − 1/p_A) (structural / D-driven component)

### Results (n = 10,000 partners)

| Phenomenon | % of partners |
|---|---|
| **C-flip** (sign of median C differs EAS vs EUR) | **32.7 %** |
| mode-class flip (most common class differs) | 18.5 % |
| any-DISJOINT flip (DISJOINT in one super-pop, not the other) | 17.7 % |

### Most informative cross-tab: among the 3,269 C-flip partners

| Sub-pattern | n | % |
|---|---|---|
| C-flip **AND** mode-class also flips | 745 | 22.8 % |
| C-flip **but** mode-class STAYS the same | 2,524 | **77.2 %** |

Lange 2013's two examples are representative of this 22.8 / 77.2 split:
- **rs2244546** (~69 kb, in window) — C-flip + class flip (CEU DISJOINT) ← 22.8 % minority pattern
- **rs9275572** (~1.3 Mb, out of window, but pattern similar) — C-flip + class STAYS PARTIAL ← 77.2 % majority pattern

### Conclusion at scale

The 2-SNP "mixed answer" was correct in spirit but understated the
prevalence of the C-driven flip. **In the MICA region, the *typical*
cross-population flip is C-driven and topology-class-silent**: signed
correlation (or equivalently, C in the M+C decomposition) sees it,
class doesn't.

- **C is the most sensitive scalar invariant** for cross-pop
  structural flip in this region.
- **Class flip is a stricter, rarer phenomenon** that requires C-flip
  *plus* the partner crossing a regime boundary (DISJOINT, STRICT,
  IDENTICAL). The latter is largely a function of the marginal
  allele-frequency mismatch *and* the magnitude |C|.
- The framework that paper 1 was building captures the rare/boundary
  cases (~20 %) rigorously, but ~80 % of the cross-pop signal in this
  region lives in the continuous C component, not in class transitions.

### Files

| File | Purpose |
|---|---|
| `data/region_500kb.vcf.gz` (+ .tbi) | 1000G phase 3, chr6:31,116,595–31,616,595, 20,019 biallelic SNVs |
| `analyze_region_MC.py` | per-partner × per-pop M/C/r/class scan |
| `render_region_figure.py` | M-C plane (panel A) + spatial scan (panel B) |
| `results/region_per_pop.parquet` | long-format per-pop data |
| `results/region_per_partner.csv` | partner-summary across pops |
| `results/figure_region_MC.{pdf,png}` | rendered figure |

## Extension 2: anchor-carrier haplotype branch decomposition (added 2026-05-08)

User's hypothesis:
> rs2596542-T does not tag a single haplotype; it tags a UNION of
> multiple haplotype branches H1, H2, ..., whose population-specific
> composition differs. Same anchor allele → different functional
> backgrounds → C-component flip across populations.

### Method

1. Subset to rs2596542-T-carrying haplotypes (n = 675 across 10 pops).
2. Build binary matrix on 7,303 SNVs with MAF ≥ 0.05 within carriers
   (±250 kb window).
3. Hierarchical clustering, Hamming distance, average linkage; cut at
   k = 2..5.
4. Per-population branch frequency table.
5. PCA visual sanity check.

### Result — hypothesis confirmed

#### Anchor carrier counts per pop (rs2596542-T = forward T = Kumar's "A risk")

| pop | n carriers |
|---|---|
| CHB 56 / JPT 69 / CHS 48 / CDX 43 / KHV 61 | EAS 277 |
| CEU 59 / TSI 89 / FIN 74 / GBR 79 / IBS 97 | EUR 398 |

#### k = 3 branch composition (frequency within rs2596542-T carriers)

| pop | H1 | H2 | H3 |
|---|---|---|---|
| **JPT** | **0.55** | 0.16 | 0.29 |
| CHB | 0.29 | 0.29 | 0.43 |
| CHS | 0.33 | 0.23 | 0.44 |
| CDX | 0.35 | 0.26 | 0.40 |
| KHV | 0.23 | 0.43 | 0.34 |
| **CEU** | **0.12** | 0.44 | 0.44 |
| TSI | 0.19 | 0.32 | 0.49 |
| **FIN** | **0.08** | 0.41 | 0.51 |
| GBR | 0.14 | 0.44 | 0.42 |
| IBS | 0.21 | 0.32 | 0.47 |

**EAS vs EUR average:**

| | H1 | H2 | H3 |
|---|---|---|---|
| EAS | 0.350 | 0.271 | 0.379 |
| EUR | 0.147 | 0.385 | 0.468 |
| Δ (EAS−EUR) | **+0.203** | −0.114 | −0.089 |

H1 frequency varies **7-fold** between JPT (0.55) and FIN (0.08) **among
carriers of the same anchor allele**.

At k = 4, an additional EUR-specific branch emerges (H3_k4 freq:
0.008 EAS / 0.154 EUR — practically EAS-absent).

#### PCA (`figure_branch_pca.png`)

PC1 = 17.8 %, PC2 = 14.1 %. Multiple visually distinct clusters; EAS
populations (blue tones) and EUR populations (green tones) populate
the clusters in unequal proportions, consistent with the branch
composition table.

### Interpretation

The hypothesis is directly supported:

1. rs2596542-T tags H1 ∪ H2 ∪ H3 — a *union* of haplotype branches,
   not a single haplotype.
2. Branch composition differs systematically across populations
   (H1 EAS-enriched; H2/H3/H4 EUR-enriched).
3. The C-component flip across populations (33 % of partners,
   region-wide) is the mathematical reflection of the continuous
   shift in branch composition: as branch frequencies change, the
   covariance between rs2596542 and each nearby SNV shifts
   smoothly — including sign reversals — without the topology class
   needing to change.
4. This explains the 77 % "C-flip without class-flip" pattern: most
   partners are correlated with multiple branches simultaneously, so
   class-level structure stays PARTIAL while the *signed* covariance
   tracks branch composition.

### Files

| File | Purpose |
|---|---|
| `cluster_anchor_haplotypes.py` | hierarchical clustering + PCA |
| `results/anchor_haplotype_branches.parquet` | haplotype-level labels (n=675) |
| `results/branch_pop_composition.csv` | per-pop branch frequencies, k=2..5 |
| `results/figure_branch_pca.{pdf,png}` | PCA scatter (pop / branch) |
| `results/figure_branch_composition.{pdf,png}` | stacked bars k=2..5 |

### What's left (deferred)

Step 4 from user's plan — annotate each branch with known functional
markers (MICA expression eQTLs from GTEx, MICA*XX coding alleles,
HLA-B tag SNPs, HCP5 markers, MICA promoter/UTR variants). This
requires external lookups (GTEx liver eQTL portal; published MICA
allele-tag SNP catalogs; HLA-B*57:01 / *58:01 etc tag tables).
Would convert the abstract H1/H2/H3 labels into "this branch is the
MICA*001 / *008 / *019 background" — interpretive payoff, not
methodological extension.

## Decision

The exploration arc is complete and the hypothesis is supported.
Results documented in this directory; no immediate publication
target. If revisited, step 4 (functional annotation overlay) is the
natural next deliverable.
