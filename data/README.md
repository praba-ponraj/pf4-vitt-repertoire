# Data

`manifest/` contains provenance metadata for the OAS inputs:

-   `bulk_oas_manifest.tsv` --- 23 bulk light-chain runs from 3 healthy
    adults (Soto et al.), 85,186,188 sequences, assigned to
    Subject-HIP1/2/3.
-   `paired_oas_manifest.tsv` --- 91 paired data units from 6
    individuals across 4 studies, 568,139 heavy--light pairs.

`derived/` contains aggregate analysis outputs used to reproduce
reported manuscript results. Run-level rows are retained in
`bulk_iglv321_results.tsv`; the manuscript reports run-level
heterogeneity in HIP2 from these rows.

Raw OAS sequence files are not included. Obtain them from OAS using the
source and accession information in the manifests; the original data
remain subject to their original sources and terms.
