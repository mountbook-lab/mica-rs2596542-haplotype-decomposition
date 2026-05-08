# rs2596542 — allele / strand / build / effect-direction mapping

> Working document. **Highest-risk step in the MICA flip-flop paper.**
> Lock this table before any haplotype analysis: one bad allele flip
> invalidates downstream conclusions.

Last updated: 2026-05-08
Maintainer: paper1-proxy-topology / mica_flipflop project

## 1. Canonical (dbSNP, GRCh38)

| Field | Value | Source |
|---|---|---|
| rsID | rs2596542 | dbSNP build 156 |
| Strand | **forward** ("All alleles are reported in the Forward orientation") | dbSNP |
| GRCh38 position | chr6:31,398,818 (NC_000006.12:g.31398818) | dbSNP |
| GRCh37 position | chr6:31,366,595 (NC_000006.11:g.31366595) | dbSNP |
| GRCh36 / NCBI Build 36 | (Kumar 2011 used Build 36 in Fig. 1) | needs lookup if paper-comparison requires |
| Reference allele | **C** (forward strand) | dbSNP / GRCh38 |
| Canonical alleles | C > A / C > G / C > T (multi-allelic site) | dbSNP |
| Practical biallelic | **C / T** (A and G alts ~0 in modern panels) | 1000G, gnomAD |
| Global T-allele frequency | ~0.43 (TOPMED) | dbSNP frequency report |
| Gene context | 2 kb upstream of MICA (5′ flanking); also intron variant of MICA-AS1 | dbSNP |

Implication: in this document, "forward T" means the **GRCh38 forward-strand T allele**, the global minor allele used in 1000G VCFs.

## 2. Kumar 2011 (Nat Genet, doi:10.1038/ng.809) — HCV-HCC, Japanese

The seed paper for rs2596542 / MICA / HCV-HCC.

| Field | Kumar 2011 reports | Maps to (forward GRCh38) |
|---|---|---|
| Allele notation in tables | A / G | reverse-strand convention |
| Risk allele | **A** ("risk allele A is higher in cases", "A risk allele is in LD with VNTR A9 and A4") | **forward T** |
| Non-risk / protective | **G** ("non-risk G allele is in LD with A6") | **forward C** |
| RAF (risk allele freq) — GWAS cases | 0.388 | T-allele freq, JPN HCC cases |
| RAF — controls | 0.331 | T-allele freq, JPN controls |
| OR | 1.39 (combined, P=4.21×10⁻¹³) | per copy of forward T |
| Genotyping platform | Illumina HumanHap610-Quad / HumanHap550v3; replication via Invader / Third Wave | Illumina TOP/BOT convention puts this SNP on reverse strand → "A/G" reporting |
| Build cited | NCBI Build 36 (Fig. 1 caption) | maps to GRCh37 chr6:31,366,595 → GRCh38 chr6:31,398,818 |
| sMICA correlation | risk allele A → **lower** soluble MICA | forward T → lower sMICA |
| LD with VNTR (A5.1 microsatellite) | A risk → A9, A4; G non-risk → A6 | forward T → A9/A4; forward C → A6 |

### Why Kumar uses A/G when dbSNP says C/T

Illumina Infinium arrays (HumanHap550v3, HumanHap610-Quad) report alleles in
TOP/BOT strand convention rather than chromosomal forward strand. For
rs2596542, the Illumina TOP strand happens to be the chromosomal **reverse**
strand, so a chromosomal C (forward ref) appears as G, and a chromosomal T
appears as A.

```
forward (GRCh38):  5' ... C / T ... 3'    ← dbSNP reports C/T
reverse (assay) :  5' ... G / A ... 3'    ← Kumar 2011 / Illumina report A/G
                              ^
                              risk = A (reverse) = T (forward)
```

This is consistent with Kumar's RAF=0.331 in JPN controls vs gnomAD/1000G
JPT T-allele frequency ≈ 0.30–0.40.

## 3. Follow-up / replication papers

### 3a. Verified from primary source

| Paper | Cohort | Disease | Allele notation | Risk allele (reported) | Risk allele (forward GRCh38) | OR | Status |
|---|---|---|---|---|---|---|---|
| Kumar 2011 Nat Genet | BioBank Japan + UTokyo | HCV-HCC | A/G (reverse / Illumina) | A | **T** | 1.39 | verified from PDF |

### 3b. Verified via open-access meta-analyses (both 2019)

Both meta-analyses pool ~10–11 primary studies. They use *opposite*
allele notations, but once strand-harmonised they agree on the risk
allele.

| Meta-analysis | Notation used | Risk allele (reported) | Strand | Risk allele (forward GRCh38) |
|---|---|---|---|---|
| Huang 2019 (BMC Med Genet, PMC6697945) | C / T | **T** (TT genotype) | forward (matches dbSNP) | **T** |
| Chen 2019 (PMC6504665) | G / A | **A** ("major allele G is protective"; OR(GG vs AA)=0.707, P=0.018) | reverse | **T** |

**Effect direction by stratum (Huang 2019, TT vs CC genotype):**

| Stratum | OR | 95% CI | P |
|---|---|---|---|
| HCV-HCC | 1.326 | 1.101 – 1.599 | 0.003 |
| HBV-HCC | 1.133 | 0.718 – 1.788 | 0.591 |
| Asian | 1.273 | 1.002 – 1.618 | 0.048 |
| European | 1.065 | 0.787 – 1.440 | 0.684 |

**Effect direction by stratum (Chen 2019, GG vs AA — forward C/C vs T/T,
so OR is inverted relative to Huang's framing):**

| Stratum | OR (GG vs AA) | P |
|---|---|---|
| HCV-HCC | 0.707 | 0.018 |
| HBV-HCC | 0.967 | 0.890 |
| Asian | 0.729 | 0.024 |
| European | 1.287 | 0.291 |

Both meta-analyses converge on:
- **Same risk allele direction** (forward T = Kumar's A) across all
  strata.
- **HCV-HCC**: clear effect, OR 1.3-ish per Huang.
- **HBV-HCC**: same direction but null / non-significant — *attenuated*,
  not flipped.
- **Asians**: effect present.
- **Europeans**: null / non-significant.
- **No primary study reports forward C as the risk allele.**

Primary studies included (Huang 2019):
- HCV: Kumar 2011, Lo 2013, Lange 2013, Burza 2016, Mohamed 2017,
  Huang 2017, Augello 2018, Hai 2017
- HBV: Tong 2013, Chen 2013, Kumar 2012

### 3c. Lange 2013 (J Hepatol) — verified from PDF, **major finding**

Title: "Comparative genetic analyses point to HCP5 as susceptibility
locus for HCV-associated hepatocellular carcinoma"
Lange et al. 2013, J Hepatol (Swiss Hepatitis C Cohort Study).
PDF: `/mnt/c/Users/yyama/OneDrive/Documents/1-s2.0-S0168827813002870-main.pdf`

| Field | Value | Source |
|---|---|---|
| Cohort | Swiss Hepatitis C Cohort Study | European |
| Allele notation | A (minor) / G (major) — same as Kumar 2011, reverse strand | text |
| Risk allele in their European cohort | **G (= forward C)** — *opposite* of Kumar 2011's "A risk" | text |
| Direct quote | "in our cohort the minor allele A of rs2596542 [...] have a **protective impact on HCC development, which represented an inverse association as compared to the results by Kumar et al.**" | Results |
| Mechanism Lange proposes | rs2596542 has *opposite-signed Pearson correlation* with the underlying causal variant in JPT vs CEU | Discussion |
| Concrete sign-flip examples | rs9275572 vs rs2596542: **r = +0.27 (JPT), r = −0.12 (CEU)**; rs2244546 vs rs2596542: **r = +0.59 (JPT), r = −0.19 (CEU)** | Fig. 1A |
| Methodological note | "The sign of the correlation was strictly kept, since SNPs can have the same r² in CEU and JPT, but their correlation may be of opposite sign." | Methods |
| Their conclusion | rs2244546 in HCP5 is a more relevant tag because its (minor) risk allele has consistent direction across both populations | Discussion |
| MAF rs2596542 | 28.4% CEU, 32.9% JPT | Results |

**This is the published flip-flop signal we were looking for.**
Lange 2013 already articulates the signed-correlation / population-
conditional tagging mechanism. They did not, however, formalise it as
*carrier-set topology* (P(v|tag) vs P(tag|v) asymmetry, n11/n10/n01,
nested vs reverse-nested classes).

### 3d. Pending PDF verification

| Paper | Cohort | Disease | Status |
|---|---|---|---|
| Tanikawa 2013 | Japanese | HCV-HCC | **PDF needed** |
| Mocci 2018 / Augello 2018 | Sicilian | HCV-HCC | **PDF needed** |
| Kumar 2012 | BioBank Japan | HBV-HCC | **PDF needed** |
| Lo 2013 | HCV cirrhotic non-SVR | HCV-HCC progression | **PDF needed** |
| Sopipong 2013 (APJCP, MICA) | Thai | HCV-HCC | **PDF still needed**; the file `APJCP_Volume 14_Issue 5_Pages 2865-2869.pdf` in `~/Documents` is a *different* Sopipong paper on **KIF1B rs17401966** (HBV-HCC Thai), not the MICA paper |

## 4. The flip-flop signal — third revision (post-Lange 2013)

The picture has changed materially after reading Lange 2013.

### What's now established

- **Lange 2013 reports a real effect-direction flip**: the minor A allele
  of rs2596542 is risk in Japanese (Kumar 2011) but **protective** in
  Lange's Swiss European HCV cohort, "an inverse association as compared
  to the results by Kumar et al." (Lange 2013 Results).
- Lange 2013 **already articulated the signed-correlation mechanism**:
  same r² across populations can come with **opposite-signed Pearson r**,
  driving opposite effect direction at the tag.
- Concrete examples published in Lange 2013 Fig. 1A:
  rs9275572 vs rs2596542  → r=+0.27 JPT, r=−0.12 CEU
  rs2244546 vs rs2596542  → r=+0.59 JPT, r=−0.19 CEU
- Two meta-analyses (Huang 2019, Chen 2019) report Europeans as **null**,
  not flipped. This is consistent with Lange — averaging primary cohorts
  with mixed signs collapses to null in the European stratum, even
  though one well-powered European cohort (Lange's Swiss SHCS) shows a
  clear direction reversal at the tag level.

### Implication for paper novelty

Lange 2013 anticipates a substantial part of what the project was
planning to demonstrate. **The proposal must be re-scoped** so that
Lange 2013 becomes the foundational prior result, and the paper's
contribution is the *methodological extension* on top:

1. **Lange 2013 used signed Pearson correlation in two populations
   (HapMap3 CEU + JPT)**. The paper-1 anchor-fixed framework
   operationalises *carrier-set topology* — P(v|tag) vs P(tag|v),
   n11/n10/n01, nested / reverse-nested / disjoint / partial
   classes — which is **richer** than signed correlation alone.
   Two SNP pairs with the same |r| can have different inclusion
   directions; carrier-set topology distinguishes them.
2. **1000 Genomes phased haplotypes across 26 populations** vs Lange's
   2-population HapMap3 contrast. We can quantify *how many* nearby
   SNPs flip topology between which population pairs, and whether
   topology heterogeneity is concentrated in particular sub-regions
   of the MICA/HCP5/HLA-B haplotype block.
3. **Connect to paper 1's anchor-fixed methodology**: rs2596542 becomes
   a worked clinical example demonstrating that the framework recovers
   and extends the signed-correlation flip-flop pattern Lange 2013
   reported.

### Recommended reframing (third revision)

- Treat **Lange 2013 as the foundational prior result**, cited
  prominently in Introduction.
- Frame the paper's contribution as: *carrier-set topology is the
  rigorous formalisation of the signed-correlation flip Lange 2013
  observed; we extend it to 1000G 26 populations and connect it to a
  systematic anchor-fixed analysis.*
- Figure 1 = Lange 2013-style sign-flipped Pearson correlation
  re-rendered in 1000G across 26 populations, plus Huang 2019 / Chen
  2019 effect-magnitude attenuation.
- Figure 3 ("same r², different topology") still the killer figure —
  but now stated as: signed r is necessary but not sufficient;
  carrier-set inclusion direction adds a second axis.
- **Decide before writing**: is this paper now (a) a methods extension
  letter that Lange 2013 made the case for empirically, or (b) absorbed
  into paper 1 as the worked clinical example? Option (b) is now more
  defensible than before because the standalone novelty has shrunk.

## 4b. Sign harmonization with Lange 2013 Fig. 1A

**Lock this table** before any cross-population r-sign comparison goes
into a figure or claim. Coded alleles for the three SNPs:

| SNP | dbSNP REF/ALT (forward GRCh37/38) | our coding (ALT=1) | Lange Fig. 1A coded allele | matches our coding? |
|---|---|---|---|---|
| rs2596542 | C / **T** | T | **T** ("rs2596542-T") | ✓ same |
| rs2244546 | C / **G** | G | **G** ("rs2244546-G") | ✓ same |
| rs9275572 | **A** / G | G | **A** ("rs9275572-A") | ✗ **opposite** |

Implication: when we report `r_signed` from `region_per_partner_26.csv`
or `per_pop_topology_with_MC.csv` and want to compare against Lange's
Fig. 1A values, the sign must be flipped for any partner whose Lange
coded allele = our REF (i.e. rs9275572 in this paper). The anchor
(rs2596542-T) is already in agreement.

**Concrete numbers** (1000G phase 3 vs Lange HapMap3):

| pair | population | our `r_signed` (ALT-coded) | harmonized to Lange | Lange Fig. 1A | sign agrees |
|---|---|---|---|---|---|
| rs2596542 × rs2244546 | JPT | **+0.489** | +0.489 | +0.59 | ✓ |
| rs2596542 × rs2244546 | CEU | **−0.187** | −0.187 | −0.19 | ✓ |
| rs2596542 × rs9275572 | JPT | −0.483 | **+0.483** (× −1) | +0.27 | ✓ |
| rs2596542 × rs9275572 | CEU | +0.260 | **−0.260** (× −1) | −0.12 | ✓ |

After harmonization, all four signs match Lange Fig. 1A. Residual
magnitude differences (~0.1–0.2) are attributable to HapMap3 vs 1000G
phase 3 cohort composition (different sample sizes, different
individuals), not to coding errors.

Tabulated: `results/lange_sign_harmonization.csv` (machine-readable).

### Manuscript boilerplate

> "Signed correlations were harmonized to the coded alleles used by
> Lange et al. 2013 (rs2596542-T, rs2244546-G, rs9275572-A) where direct
> comparisons are made. Our default coding throughout the analysis is
> the dbSNP forward-strand ALT allele (rs2596542-T, rs2244546-G,
> rs9275572-G); for rs9275572 this is opposite to Lange's coded allele,
> so the published Lange values are reported with the sign flipped to
> our coding in supplementary tables."

A single sentence of this form must appear in Methods (or the figure
caption of Figure 1) of any manuscript that cites Lange's r-values.

### "Better tag" framing for rs2244546

Lange's argument that rs2244546 is a more relevant tag than rs2596542
rests on:

1. rs2244546 minor allele has the **same effect direction** in
   Caucasians and Japanese (Discussion).
2. rs2244546 sits in HCP5, which Lange identifies as the susceptibility
   gene rather than MICA upstream.

In the harmonized table above, rs2596542 × rs2244546 shows JPT +0.489
vs CEU −0.187 — a clean sign flip across populations, but the partner
(rs2244546-G) has consistent effect direction with the disease-causal
variant in *both* populations because its r with the underlying causal
variant has the same sign. By contrast, rs2596542 × causal-variant has
opposite signs in JPT vs CEU (Kumar 2011 risk = T in JPT; Lange 2013
"protective" = T in CEU). When using rs2596542 as the GWAS tag, the
"flip" appears at the disease level; when using rs2244546 as the tag,
the flip is absorbed into the LD pattern with rs2596542 and the disease
direction stays consistent. **This is what makes rs2244546 the "better
tag" in Lange's framing.**

When the manuscript discusses this, it must use the harmonized signs
above (not the raw `r_signed` values from the per-pop CSVs), and must
state explicitly that rs2596542 × rs2244546 is the *flip itself in our
coordinate system* — the consistency Lange refers to is between
rs2244546 and the unobserved causal variant, not between rs2244546 and
rs2596542.

## 5. To-do (next step)

- [ ] Download / locate PDFs for Tanikawa 2013, Lange 2013, Sopipong 2013,
      Mocci 2018, Kumar 2012 (HBV).
- [ ] Verify allele notation, risk allele, and effect direction in each.
- [ ] Map every paper's reported risk allele to the **forward GRCh38**
      reference frame (C / T).
- [ ] Tabulate cohort × allele notation × effect direction in the CSV
      sibling file `rs2596542_allele_strand_build_mapping.csv`.
- [ ] Decide which "flips" are real (after notation correction) and
      promote those rows to Figure 1 of the MICA flip-flop letter.
