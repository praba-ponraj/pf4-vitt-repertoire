# Healthy human B-cell repertoire reference frequencies for vaccine-induced immune thrombocytopenia and thrombosis--associated antibody sequence features

**Archived release:** [Zenodo concept DOI:
10.5281/zenodo.23087012](https://doi.org/10.5281/zenodo.23087012)
(resolves to the current version; v1.1.0 adds the quality-control and
sensitivity analyses reported in the manuscript)

## Repository purpose

This repository contains analysis code, provenance metadata, and derived
aggregate outputs supporting the manuscript *Healthy human B-cell
repertoire reference frequencies for vaccine-induced immune
thrombocytopenia and thrombosis--associated antibody sequence features*.

The journal does not permit Supporting Information for Brief Reports.
This repository is therefore the sole location of the extended
methodological detail, and is organized so that **every quantitative
statement in the manuscript can be traced to the script and output field
that produced it** (see *Manuscript traceability* below).

## Version history

  -----------------------------------------------------------------------
  Version                             Contents
  ----------------------------------- -----------------------------------
  1.0.0                               Bulk and paired IGLV3-21 analyses;
                                      OAS provenance manifests; derived
                                      aggregate outputs.

  1.1.0                               Adds CDR-L3 distribution and
                                      abundance-removal sensitivity
                                      analyses, run-level quality
                                      control, and the sequence-observed
                                      DDS sensitivity analysis. Separates
                                      analyses examined but not reported.
                                      Corrects the K31E positional
                                      description (see *Feature
                                      definitions*).
  -----------------------------------------------------------------------

## Data source

The primary repertoire data are publicly available through Observed
Antibody Space (OAS). Raw OAS sequence files are not redistributed here.
Source and accession provenance are recorded in the manifest files under
`data/manifest/`; users should obtain the corresponding datasets from
OAS. Original OAS and source data remain subject to their original
sources and terms.

## Repository structure

``` text
.
├── README.md
├── CITATION.cff
├── LICENSE
├── scripts/
│   ├── analyze_bulk_iglv321.py        # bulk light-chain analysis (v1.0.0)
│   ├── analyze_paired_iglv321.py      # paired analysis (v1.0.0)
│   ├── oas_manifest.py                # provenance manifest utility (v1.0.0)
│   ├── k31e_clonality.py              # CDR-L3 distribution (v1.1.0)
│   ├── k31e_allele_runs.py            # abundance-removal + run-level QC (v1.1.0)
│   └── paired_dds_sensitivity.py      # sequence-observed DDS subset (v1.1.0)
├── data/
│   ├── README.md
│   ├── manifest/
│   │   ├── bulk_oas_manifest.tsv
│   │   └── paired_oas_manifest.tsv
│   └── derived/
│       ├── bulk_iglv321_results.tsv
│       └── paired_iglv321_results.tsv
└── results/
    ├── reported/                      # outputs supporting reported results
    │   ├── clon_HIP1.txt
    │   ├── clon_HIP2.txt
    │   ├── clon_HIP3.txt
    │   ├── HIP1_K31E_cdr3.tsv
    │   ├── HIP2_K31E_cdr3.tsv
    │   ├── HIP3_K31E_cdr3.tsv
    │   ├── qc_HIP1_all.txt
    │   ├── qc_HIP2_all9.txt
    │   ├── qc_HIP3_all.txt
    │   └── sens_table1c.txt
    └── exploratory_not_reported/
        ├── README.md
        ├── allele_association/
        └── hcdr3_motif_search/
```

## Manuscript traceability

Every number in the manuscript maps to one script and one output field.

| Manuscript location | Claim | Script | Output file | Field / row |
|---|---|---|---|---|
| Abstract, Results | K31E in 43,585 of 1,455,792 IGLV3-21 sequences (2.99%) | `analyze_bulk_iglv321.py` | `data/derived/bulk_iglv321_results.tsv` | `TOTAL` row, `k31e_all_iglv321` / `iglv321` |
| Results, Table 1 | HIP1 18,727/545,394 (3.43%); HIP2 6,564/328,217 (2.00%); HIP3 18,294/582,181 (3.14%) | `analyze_bulk_iglv321.py` | `data/derived/bulk_iglv321_results.tsv` | subject-level rows |
| Methods | 773 of 1,455,792 (0.053%) could not be anchored | `analyze_bulk_iglv321.py` | `data/derived/bulk_iglv321_results.tsv` | `TOTAL` row, `k31_unmapped_iglv321` |
| Table 1 | \*02/\*03 assignment proportions per donor | `analyze_bulk_iglv321.py` | `data/derived/bulk_iglv321_results.tsv` | `allele_02_03_pct`, subject rows |
| Results | K31E across 1,093 to 3,778 distinct CDR-L3 amino-acid sequences per individual | `k31e_clonality.py` | `results/reported/clon_HIP1.txt`, `clon_HIP2.txt`, `clon_HIP3.txt` | `unique CDR3`, K31E column |
| Results | Per-CDR-L3 K31E counts underlying CDR-L3 diversity and abundance-removal checks | `k31e_clonality.py` | `results/reported/HIP1_K31E_cdr3.tsv`, `HIP2_K31E_cdr3.tsv`, `HIP3_K31E_cdr3.tsv` | CDR-L3 sequence and K31E count rows |
| Results | K31E remained at 1.12% to 2.21% after removal of the 25 most abundant CDR-L3 sequences | `k31e_allele_runs.py` | `results/reported/qc_HIP1_all.txt`, `qc_HIP2_all9.txt`, `qc_HIP3_all.txt` | section (A), `top 25` row, `prevalence` |
| Results | HIP2 run-level K31E heterogeneity (0.16%–2.17%) | `analyze_bulk_iglv321.py` | `data/derived/bulk_iglv321_results.tsv` | run-level rows, `k31e_pct` |
| Table 2, column 1 | 9,394 proxy pairs; K31E 133 (1.42%); IGHJ3 1,613 (17.17%); acidic HCDR3 2,556 (27.21%); VITT IGHV 1,839 (19.58%); VITT IGHV + IGHJ3 309 (3.29%); five-feature 0 | `analyze_paired_iglv321.py` or `paired_dds_sensitivity.py` | `data/derived/paired_iglv321_results.tsv`; `results/reported/sens_table1c.txt` | `TOTAL` row; `full proxy` column |
| Methods, Table 2 col. 2 | 7,424 of 9,394 (79.0%) with DDS directly observed; K31E 52 (0.70%) | `paired_dds_sensitivity.py` | `results/reported/sens_table1c.txt` | `DDS observed` column |
| Results, Table 2 col. 3 | 1,970 (21.0%) qualifying by \*02/\*03 assignment alone; K31E 81 (4.11%) | `paired_dds_sensitivity.py` | `results/reported/sens_table1c.txt` | `allele-only` column |

Provenance chain: OAS source/accession → manifest → analysis script →
derived output → manuscript table.

## Input data

Use the provenance information in the manifests to obtain the
corresponding OAS `.csv.gz` files. Raw inputs should remain outside this
repository.

For the bulk analysis, organize the OAS files into three
subject-specific subdirectories:

``` text
INPUT_DIR/
├── Subject-HIP1/
│   └── *.csv.gz
├── Subject-HIP2/
│   └── *.csv.gz
└── Subject-HIP3/
    └── *.csv.gz
```

The assignment of OAS accessions to HIP1, HIP2, and HIP3 is provided in
`data/manifest/bulk_oas_manifest.tsv`.

For the paired analysis, place the OAS `.csv.gz` files directly in a
single input directory. The required 91 data units are listed in
`data/manifest/paired_oas_manifest.tsv` and total 568,139 paired
sequences.

To generate a manifest from one or more directories (searched
recursively), write standard output to a file:

``` bash
python3 scripts/oas_manifest.py DIR1 DIR2 DIR3 > oas_manifest.tsv
```

Warnings and diagnostics are written to standard error and therefore do
not corrupt the TSV stream.

## Feature definitions

**K31E** was defined as glutamate at the sequence position aligned to
the IGLV3-21 germline lysine. The position is located **per sequence**
by searching `v_germline_alignment_aa` for the germline context `NIGSKS`
and taking offset +4, with the germline residue at that offset verified
to be lysine before the call is made. No fixed query index or IMGT
coordinate is used. Because the anchor is sought in the germline
alignment rather than in the query, the definition is unaffected by 5′
truncation or by somatic mutation elsewhere in the sequence. Sequences
in which the anchor cannot be located are counted separately and
excluded from the numerator.

*Note (v1.1.0):* v1.0.0 described this position as "IMGT position 37".
That coordinate was not independently verified and has been removed. The
operational definition above is unchanged and is what the code
implements.

For paired IGLV3-21 sequences, the **DDS(D) proxy** was defined as
IGLV3-21 assigned to `*02` or `*03`, **or** containing DDS in the OAS
`cdr2_aa_light` field. The `*02` and `*03` alleles encode DDS in CDR-L2
followed by an adjacent D residue across the CDR-L2/framework boundary.
Because this is a disjunction, `paired_dds_sensitivity.py` reports the
proxy population split into pairs with DDS directly observed and pairs
qualifying by allele assignment alone.

An **acidic HCDR3** contained at least three aspartate or glutamate
residues (`D + E >= 3`). The **VITT-reported IGHV genes** were
`IGHV3-7`, `IGHV5-51`, `IGHV2-26`, `IGHV3-30`, and `IGHV7-4-1`.

Allele assignments reported here are **descriptive repertoire
annotations derived from OAS, not genomic genotypes.** No genotype is
inferred and no allele association is claimed.

## Requirements

Python 3.

  Script                        Dependencies
  ----------------------------- -----------------------
  `analyze_bulk_iglv321.py`     standard library only
  `analyze_paired_iglv321.py`   standard library only
  `oas_manifest.py`             standard library only
  `paired_dds_sensitivity.py`   standard library only
  `k31e_clonality.py`           `pandas`
  `k31e_allele_runs.py`         `pandas`

## Reproduction commands

Bulk analysis (v1.0.0):

``` bash
python3 scripts/analyze_bulk_iglv321.py INPUT_DIR \
    --output bulk_iglv321_results.tsv
```

Paired analysis (v1.0.0):

``` bash
python3 scripts/analyze_paired_iglv321.py INPUT_DIR \
    --output paired_iglv321_results.tsv
```

CDR-L3 distribution of K31E (v1.1.0), one donor per invocation:

``` bash
python3 scripts/k31e_clonality.py Subject-HIP1/*.csv.gz --label HIP1
python3 scripts/k31e_clonality.py Subject-HIP2/*.csv.gz --label HIP2
python3 scripts/k31e_clonality.py Subject-HIP3/*.csv.gz --label HIP3
```

Abundance-removal sensitivity and run-level quality control (v1.1.0):

``` bash
python3 scripts/k31e_allele_runs.py Subject-HIP1/*.csv.gz --label HIP1
python3 scripts/k31e_allele_runs.py Subject-HIP2/*.csv.gz --label HIP2_all9
python3 scripts/k31e_allele_runs.py Subject-HIP3/*.csv.gz --label HIP3
```

Sequence-observed DDS sensitivity analysis (v1.1.0):

``` bash
python3 scripts/paired_dds_sensitivity.py PAIRED_INPUT_DIR/*.csv.gz
```

## Expected reproducibility results

### Bulk light-chain analysis

| Group | K31E / IGLV3-21 | Frequency |
|---|---:|---:|
| HIP1 | 18,727 / 545,394 | 3.434% |
| HIP2 | 6,564 / 328,217 | 2.000% |
| HIP3 | 18,294 / 582,181 | 3.142% |
| Total | 43,585 / 1,455,792 | 2.994% |

The bulk input contained 85,186,188 sequences; 773 IGLV3-21 sequences
(0.053%) had an unanchored K31 position.

### CDR-L3 distribution of K31E

The distinct-CDR-L3 count for each donor is reported independently by
`k31e_clonality.py` and by section (A) of `k31e_allele_runs.py`; the two
agree exactly (e.g. HIP2: 1,093 in both).

| Donor | K31E sequences | Distinct CDR-L3 sequences | Most abundant CDR-L3, share |
|---|---:|---:|---:|
| HIP1 | 18,727 | 3,778 | 6.37% |
| HIP2 (all 9 runs) | 6,564 | 1,093 | 13.54% |
| HIP3 | 18,294 | 2,762 | 8.42% |

### Abundance-removal sensitivity

K31E frequency among anchored IGLV3-21 sequences after removing the N
most abundant K31E-associated CDR-L3 amino-acid sequences:

| Removed | HIP1 | HIP2 (all 9 runs) | HIP3 |
|---|---:|---:|---:|
| none | 3.4365% | 2.0005% | 3.1435% |
| top 1 | 3.2176% | 1.7296% | 2.8789% |
| top 5 | 2.8517% | 1.4401% | 2.4498% |
| top 10 | 2.6454% | 1.3300% | 2.1860% |
| top 25 | 2.2101% | 1.1240% | 1.8017% |

### Run-level heterogeneity (HIP2)

Among HIP2 runs with at least 1,000 IGLV3-21 sequences, K31E frequency ranged
from 0.16% to 2.17%. Three runs with the highest IGLV3-21 yield (~2.3% of
sequences) gave 2.06–2.17%; three with lower yield (~0.29%) gave 0.16–0.25%.
Three further runs yielded 6, 11 and 12 IGLV3-21 sequences and contribute
negligibly. This heterogeneity was identified during analysis and is reported
as a limitation; no run was excluded.

### Paired analysis and sequence-observed DDS sensitivity

| Measure | Full DDS(D) proxy | DDS observed | \*02/\*03 only |
| --- | ---: | ---: | ---: |
| Pairs | 9,394 | 7,424 (79.03%) | 1,970 (20.97%) |
| K31E | 133 (1.42%) | 52 (0.70%) | 81 (4.11%) |
| IGHJ3 | 1,613 (17.17%) | 1,300 (17.51%) | 313 (15.89%) |
| Acidic HCDR3 (≥3 D/E) | 2,556 (27.21%) | 2,033 (27.38%) | 523 (26.55%) |
| VITT-reported IGHV | 1,839 (19.58%) | 1,464 (19.72%) | 375 (19.04%) |
| VITT IGHV + IGHJ3 | 309 (3.29%) | 237 (3.19%) | 72 (3.65%) |
| Complete five-feature combination | 0 | 0 | 0 |

Total paired input: 568,139 pairs from 91 OAS data units; 16,464
IGLV3-21 pairs; 0 unanchored K31 positions.

## Analyses examined but not reported

`results/exploratory_not_reported/` contains analyses performed during
manuscript development that are **not** used as evidence for any
quantitative result in the manuscript. See that directory's README for
what was examined and why it was not reported.

## License

Repository code is provided under the MIT License. This license does not
apply to third-party OAS or source data, which remain subject to their
original sources and terms.
