# VITT-Associated Antibody Features in Healthy Human B-Cell Repertoires
**Archived release:** [Zenodo DOI: 10.5281/zenodo.23087012](https://doi.org/10.5281/zenodo.23087012)

## Repository purpose

This repository contains analysis code, provenance metadata, and derived aggregate outputs supporting the manuscript *VITT-Associated Antibody Features in Healthy Human B-Cell Repertoires*.

## Data source

The primary repertoire data are publicly available through Observed Antibody Space (OAS). Raw OAS sequence files are not redistributed here. Source and accession provenance are recorded in the manifest files under `data/manifest/`; users should obtain the corresponding datasets from OAS. Original OAS and source data remain subject to their original sources and terms.

## Repository structure

- `scripts/`: validated bulk and paired analysis scripts, plus the OAS manifest utility.
- `data/manifest/`: source and accession provenance for the OAS inputs.
- `data/derived/`: derived aggregate outputs supporting the reported manuscript results.

## Input data

Use the provenance information in the manifests to obtain the corresponding OAS `.csv.gz` files. Raw inputs should remain outside this repository.

For the bulk analysis, organize the OAS files into three subject-specific subdirectories:

```text
INPUT_DIR/
├── Subject-HIP1/
│   └── *.csv.gz
├── Subject-HIP2/
│   └── *.csv.gz
└── Subject-HIP3/
    └── *.csv.gz
```

The assignment of OAS accessions to HIP1, HIP2, and HIP3 is provided in `data/manifest/bulk_oas_manifest.tsv`.

For the paired analysis, place the OAS `.csv.gz` files directly in a single input directory. The required accessions are provided in `data/manifest/paired_oas_manifest.tsv`.

To generate a manifest from one or more directories (searched recursively), write standard output to a file:

```bash
python3 scripts/oas_manifest.py DIR1 DIR2 DIR3 > oas_manifest.tsv
```

Warnings and diagnostics are written to standard error and therefore do not corrupt the TSV stream. The repository includes separate provenance manifests for the bulk (`data/manifest/bulk_oas_manifest.tsv`) and paired (`data/manifest/paired_oas_manifest.tsv`) datasets.

## Feature definitions

K31E was defined as glutamate at the sequence position aligned to the IGLV3-21 germline lysine corresponding to IMGT position 37 (K31 in Wang et al.). This alignment-anchored definition is used for both the bulk and paired analyses.

For paired IGLV3-21 sequences, the DDS(D) proxy was defined as IGLV3-21 assigned to `*02` or `*03`, or containing DDS in CDR-L2. The `*02` and `*03` alleles encode DDS in CDR-L2 followed by an adjacent D residue across the CDR-L2/framework boundary.

An acidic HCDR3 contained at least three aspartate or glutamate residues (`D + E >= 3`). The VITT-reported IGHV genes were `IGHV3-7`, `IGHV5-51`, `IGHV2-26`, `IGHV3-30`, and `IGHV7-4-1`.

## Requirements

The public scripts use only the Python standard library. A Python 3 interpreter is required; no third-party Python packages are needed.

## Reproduction commands

Bulk analysis:

```bash
python3 scripts/analyze_bulk_iglv321.py INPUT_DIR \
    --output bulk_iglv321_results.tsv
```

Paired analysis:

```bash
python3 scripts/analyze_paired_iglv321.py INPUT_DIR \
    --output paired_iglv321_results.tsv
```

## Expected reproducibility results

Bulk results:

| Group | K31E / IGLV3-21 | Frequency |
| --- | ---: | ---: |
| HIP1 | 18,727 / 545,394 | 3.434% |
| HIP2 | 6,564 / 328,217 | 2.000% |
| HIP3 | 18,294 / 582,181 | 3.142% |
| Total | 43,585 / 1,455,792 | 2.994% |

The bulk input contained 85,186,188 sequences; 773 IGLV3-21 sequences had an unmapped K31 position.

Paired results:

| Measure | Result |
| --- | ---: |
| Total paired BCRs | 568,139 |
| IGLV3-21 | 16,464 |
| K31E among all IGLV3-21 | 278 (1.689%) |
| DDS(D) proxy | 9,394 |
| K31E among proxy | 133 (1.42%) |
| IGHJ3 | 1,613 (17.2%) |
| Acidic HCDR3 | 2,556 (27.2%) |
| VITT-reported IGHV | 1,839 (19.6%) |
| VITT-reported IGHV + IGHJ3 | 309 (3.3%) |
| Complete five-feature combination | 0 |
| K31 unmapped among paired IGLV3-21 | 0 |

The complete combination evaluates the DDS(D) proxy, K31E, VITT-reported IGHV, IGHJ3, and acidic HCDR3.

## Manuscript table mapping

- `data/derived/bulk_iglv321_results.tsv` supports and reproduces manuscript Table 1A.
- `data/derived/paired_iglv321_results.tsv` supports and reproduces manuscript Table 1C.

## Provenance chain

OAS source/accession → manifest → analysis script → derived TSV → manuscript table.

## License

Repository code is provided under the MIT License. This license does not apply to third-party OAS or source data, which remain subject to their original sources and terms.
