# Data Sources and Redistribution

This workflow reuses public cancer research resources. The repository
publishes code, configuration, provenance records, and small derived summary
tables. It does not redistribute raw expression counts, clinical source files,
or DepMap matrices. Obtain those from their providers under the providers'
terms.

## Resources Used

| Resource | What is used | Access tier | Version or snapshot | How it enters the workflow |
| --- | --- | --- | --- | --- |
| NCI Genomic Data Commons (GDC), project TCGA-BRCA | STAR gene-count files for primary tumors and solid-tissue normals (1,224 files) | Open access | Frozen local download; file-to-sample metadata retrieved 2026-08-20T15:25:55Z | Read from `GDCdata/` (not committed). One metadata-only lookup mapped file UUIDs to samples: `analysis_v2/provenance/gdc_expression_metadata.tsv` and `.query.json`. |
| GDC TCGA-BRCA clinical supplement | `nationwidechildrens.org_clinical_patient_brca.txt` (receptor status, age, race, tissue source site) | Open access | Download manifest with MD5 values: `MANIFEST.txt` | Read from `GDCdata/` (not committed). |
| DepMap | `CRISPRGeneEffect.csv`, `CRISPRInferredCommonEssentials.csv`, `Gene.csv`, `Model.csv`, `ScreenGeneDependency.csv`, `ScreenGeneEffect.csv` | Public release | Unresolved local snapshot (warning W001); files are hashed in `analysis_v2/provenance/input_checksums.sha256` | Read from `DepMapCrispr/` (not committed). |
| Gene Ontology through Bioconductor | `org.Hs.eg.db` and `GO.db` | Public | 3.20.0 | Installed R packages; no online query. |
| MSigDB KEGG legacy collection through `msigdbr` | C2:CP:KEGG_LEGACY gene sets | Public, subject to MSigDB and KEGG terms | MSigDB 2024.1.Hs | Installed R package; no online query. |
| Curated tractability evidence | None supplied | Not applicable | Not applicable | Empty schema in `analysis_v2/config/tractability_evidence.csv`; tractability is reported as not assessed. |

No controlled-access (dbGaP-authorized) data were used.

## What This Repository Redistributes

- GDC file metadata linking open-access file UUIDs to TCGA sample barcodes.
- Checksums of the frozen inputs, so a user who downloads the same files can
  confirm they match.
- Legacy summary tables in `results/` and `DepMapCrispr/depmap_integration_outputs/`
  derived from TCGA and DepMap. These are retained as historical evidence of
  the audited exploratory workflow and are not the corrected results.

## What It Does Not Redistribute

- Raw STAR count files or clinical source files from GDC.
- DepMap source matrices.
- Gene-set collection files.
- Generated corrected outputs under `analysis_v2/results/` and
  `analysis_v2/plots/`, which are excluded by `.gitignore`.

## Obtaining the Inputs

1. Download the TCGA-BRCA open-access STAR count files and the clinical
   supplement listed in `MANIFEST.txt` from the GDC, for example with
   `TCGAbiolinks`, and place them under `GDCdata/TCGA-BRCA/` using the paths in
   `analysis_v2/config/analysis.yml`.
2. Download the six DepMap files listed above from the DepMap portal into
   `DepMapCrispr/`.
3. From the repository root, confirm the downloaded bytes match the frozen
   snapshot before running anything else:
   `sha256sum -c analysis_v2/provenance/input_checksums.sha256`.
   Do this first, because stage 00 currently rewrites that file rather than
   comparing against it (finding O1 in `docs/AUDIT_REMEDIATION.md`).
4. Run `Rscript analysis_v2/scripts/00_verify_inputs.R`. It fails if required
   inputs are missing or the local expression files do not match the frozen GDC
   metadata.

GDC files can be reprocessed or retired by the provider, and the current
DepMap release will differ from the unresolved local snapshot. A byte-identical
reconstruction of the inputs is therefore not guaranteed.

## Attribution

Results derived from these resources should cite the TCGA Research Network and
the NCI GDC, the DepMap project, the Gene Ontology Consortium, and MSigDB and
KEGG, following each provider's current citation guidance.
