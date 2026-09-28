# Corrected TNBC Analysis Design

**Status:** Approved in discussion; awaiting review of this written specification  
**Date:** 2026-08-19

## Purpose

Create a reproducible, separately versioned correction of the exploratory TNBC target-discovery workflow. The corrected analysis must preserve the current TCGA/GDC and DepMap snapshots, repair cohort and gene-identity errors, strengthen the statistical design, and keep scientific evidence streams separate. The legacy script and its outputs remain unchanged as historical exploratory work.

## Selected Approach

The corrected workflow will be implemented in parallel under `analysis_v2/`. This was selected over an in-place rewrite, which would obscure the provenance of existing results, and a minimal patch, which would leave major statistical and target-ranking weaknesses unresolved.

```text
analysis_v2/
├── R/                 # Reusable analysis and validation functions
├── scripts/           # Numbered executable stages
├── config/            # Explicit thresholds and settings
├── provenance/        # Checksums, manifests, versions, and run metadata
├── tests/testthat/     # Unit and regression tests
├── results/           # Corrected, versioned result tables
└── plots/             # Corrected, versioned figures
```

The workflow will not overwrite legacy `results/`, `plots/`, or `tnbc_analysis_claude.r` artifacts.

## Frozen Inputs and Provenance

The analysis will use the currently downloaded TCGA/GDC and DepMap files without replacing or expanding those data snapshots. Because the downloaded expression files do not contain sample barcodes and no expression query manifest was retained, one metadata-only GDC API lookup is permitted to map the existing file UUIDs to TCGA sample identifiers. The response will be frozen with its retrieval time and checksum and reused offline thereafter. No expression or clinical data will be downloaded by this lookup. Before computation, the workflow will inventory required inputs, record file checksums and available source identifiers, and stop if a required file is missing or changes unexpectedly.

Every completed run will record its command, configuration, input and output manifests, package versions, `sessionInfo()`, warnings, and cohort flow. A source release that cannot be established from the local files will be labeled unresolved rather than inferred.

## Cohort Construction

The primary contrast is eligible primary TNBC tumors versus all solid-tissue normal breast samples.

- Sample type will come from TCGA metadata, with barcode-derived values used only as validated supporting information.
- ER, PR, and HER2 evidence will determine TNBC status for primary tumors only.
- Normal-tissue eligibility will be independent of the associated patient's receptor fields.
- One biological sample per patient will be chosen through a documented deterministic rule.
- Ambiguous receptor classifications will not silently enter the TNBC group.
- A sample manifest will state patient, sample, group, receptor evidence, pairing status, and every inclusion or exclusion reason.

The sample manifest is the authoritative input to downstream contrasts.

## Differential Expression

Version-stripped Ensembl identifiers will remain the primary gene keys. Gene symbols will be annotations; duplicate symbols will not be altered with synthetic suffixes.

Low-expression filtering will be calculated from the samples participating in each contrast. The primary DESeq2 model will estimate TNBC versus normal effects. Batch or clinical covariates will be added only after documented completeness, variation, and non-confounding checks. Effect-size ranking will use shrunken log-fold changes, while inferential results will retain adjusted p-values and the unshrunk model estimates needed for audit.

Two sensitivity analyses will be reported separately:

1. A paired TNBC-normal analysis with patient represented in the design.
2. TNBC versus other primary breast tumors to assess TNBC specificity.

The report will show effect direction, uncertainty, and concordance across analyses. Sensitivity results will not be collapsed into the primary model or a composite score.

## Enrichment Analysis

The corrected workflow will perform both over-representation analysis on declared significant sets and ranked-list enrichment using the DE statistics. For each database, the tested and successfully mapped genes will define the eligible background. Outputs will record the tested universe, mapping success and losses, database/version information available locally, and multiple-testing adjustment.

Upregulated and downregulated results will remain separate. Missing identifiers and empty result sets will produce explicit diagnostic outputs rather than misleading blank figures.

## DepMap Integration and Candidate Evidence

MDA-MB-231 will be retained as an index model, not treated as a proxy for all TNBC. DepMap evidence will summarize a predefined, provenance-recorded TNBC cell-line panel and compare it with non-TNBC breast models when the local metadata supports those assignments.

Candidate tables will retain distinct evidence columns for:

- primary expression effect and FDR;
- paired-analysis direction and stability;
- TNBC specificity relative to other breast tumors;
- dependency across TNBC models;
- dependency difference from non-TNBC breast models;
- common-essential status; and
- curated tractability evidence with an explicit source and version.

Regex-based druggability labels and arithmetic combinations of incompatible statistics will be removed. Transparent evidence tiers may be assigned from documented rules, but raw evidence will remain visible and no tier will be described as validated therapeutic efficacy.

## Execution and Failure Behavior

The pipeline will expose restartable stages:

```text
00_verify_inputs.R
01_build_cohort.R
02_run_differential_expression.R
03_run_enrichment.R
04_integrate_depmap.R
05_build_report.R
```

Each stage will validate its input schema and upstream manifest before writing outputs. It will fail with a specific message for checksum mismatches, missing columns, non-unique primary keys, invalid cohort states, unsupported model metadata, or inconsistent result counts. Partial outputs will not be presented as a completed run.

## Validation Strategy

Fixture-based tests will cover receptor normalization and classification, normal-sample independence, duplicate-patient selection, gene-key preservation, contrast-specific filtering, and enrichment-universe construction. Regression tests will specifically prevent receptor fields from relabeling normal tissue.

Integration checks will verify required schemas, unique identifiers, nonempty expected outputs, internally consistent sample/gene counts, and recorded checksums. R scripts will receive parse checks before fixture tests and the full frozen-data run.

## Reporting and Scientific Boundaries

The final report will distinguish primary findings, sensitivity-supported findings, and exploratory candidates. It will include sample flow, exclusions, mapping losses, threshold definitions, warnings, and provenance links. Limitations will include unmatched normal tissue, incomplete clinical covariates, any unresolved data release, DepMap model annotation limits, and cell-line-to-tumor translation. The corrected workflow supports target nomination, not target validation or therapeutic claims.
