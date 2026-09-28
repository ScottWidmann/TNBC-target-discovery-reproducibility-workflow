# Corrected TNBC Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a separately versioned, reproducible TNBC target-nomination workflow from the frozen local TCGA/GDC and DepMap snapshots.

**Architecture:** Implement focused R modules under `analysis_v2/R/`, called by six restartable stage scripts. Freeze the one permitted GDC metadata lookup, make a sample manifest the cohort authority, retain Ensembl IDs through DE and enrichment, and combine expression and DepMap results as separate evidence fields rather than a composite score.

**Tech Stack:** R 4.4.2; data.table, dplyr, yaml, digest, jsonlite, DESeq2, apeglm, clusterProfiler, org.Hs.eg.db, GO.db, msigdbr, fgsea, ggplot2, ggrepel, rmarkdown, testthat, TCGAbiolinks.

**Spec:** `docs/superpowers/specs/2026-08-19-corrected-tnbc-analysis-design.md`

## Global Constraints

- Leave `tnbc_analysis_claude.r`, legacy `results/`, legacy `plots/`, and the downloaded data files unchanged.
- Use the 1,224 expression files already under `GDCdata/`; do not download or replace expression or clinical data.
- Permit exactly one metadata-only GDC query to map existing file UUIDs, then reuse the checksum-frozen response offline.
- Use all eligible primary TNBC tumors versus all solid-tissue normals as the primary comparison.
- Select one sample per patient and analysis group with a documented deterministic rule.
- Keep version-stripped Ensembl IDs as primary gene keys; gene symbols are non-unique annotations.
- Keep primary, paired, specificity, dependency, common-essential, and tractability evidence in separate columns.
- Do not use gene-name druggability regexes or an arithmetic priority score.
- Use two-space indentation, UTF-8, relative project paths, and explicit package qualification where ambiguity is possible.
- Do not silently substitute a smaller dataset, omit a failed analysis, or present partial outputs as a completed run.

## File Map

- `analysis_v2/config/analysis.yml`: paths, thresholds, seeds, and declared evidence rules.
- `analysis_v2/config/tractability_evidence.csv`: explicit empty-by-default schema for future versioned evidence.
- `analysis_v2/R/io_contracts.R`: configuration, path, checksum, schema, and GDC metadata functions.
- `analysis_v2/R/expression_io.R`: STAR-count parsing and matrix assembly.
- `analysis_v2/R/cohort.R`: clinical parsing, receptor classification, sample inclusion, and deduplication.
- `analysis_v2/R/differential_expression.R`: filtering, covariate gates, DESeq2 models, and shrinkage.
- `analysis_v2/R/enrichment.R`: identifier mapping, tested universes, ORA, and ranked enrichment.
- `analysis_v2/R/depmap.R`: model grouping, matrix extraction, dependency summaries, and common essentials.
- `analysis_v2/R/evidence.R`: evidence joins and transparent evidence classes.
- `analysis_v2/R/reporting.R`: output checks, plots, provenance summaries, and report inputs.
- `analysis_v2/scripts/00_verify_inputs.R` through `05_build_report.R`: restartable executable stages.
- `analysis_v2/corrected_analysis_report.Rmd`: final report template.
- `analysis_v2/tests/testthat/`: fixture-based unit and regression tests mirroring the modules.

---

### Task 1: Configuration, input contracts, and frozen GDC metadata

**Files:**
- Create: `analysis_v2/config/analysis.yml`
- Create: `analysis_v2/config/tractability_evidence.csv`
- Create: `analysis_v2/R/io_contracts.R`
- Create: `analysis_v2/scripts/00_verify_inputs.R`
- Create: `analysis_v2/provenance/README.md`
- Create: `analysis_v2/tests/testthat.R`
- Create: `analysis_v2/tests/testthat/test-io-contracts.R`
- Modify: `.gitignore`

**Interfaces:**
- Produces: `read_analysis_config(path) -> list`
- Produces: `discover_expression_files(root) -> data.frame(file_id, file_name, path)`
- Produces: `fetch_gdc_file_metadata(local_files, query_fn) -> data.frame`
- Produces: `validate_expression_metadata(local_files, metadata) -> invisible(TRUE)`
- Produces: `write_checksum_manifest(paths, output_file) -> data.frame`
- Produces: `require_columns(x, required, object_name) -> invisible(TRUE)`

- [ ] **Step 1: Write failing contract tests**

```r
test_that("expression discovery derives file UUIDs from frozen paths", {
  root <- tempfile("counts-")
  dir.create(file.path(root, "file-uuid"), recursive = TRUE)
  file.create(file.path(root, "file-uuid", "counts.tsv"))
  found <- discover_expression_files(root)
  expect_identical(found$file_id, "file-uuid")
  expect_identical(found$file_name, "counts.tsv")
})

test_that("metadata validation rejects missing and duplicate file IDs", {
  local <- data.frame(file_id = c("a", "b"), file_name = c("a.tsv", "b.tsv"))
  expect_error(validate_expression_metadata(local, data.frame(file_id = "a")), "missing")
  duplicate <- data.frame(file_id = c("a", "a", "b"))
  expect_error(validate_expression_metadata(local, duplicate), "duplicate")
})
```

- [ ] **Step 2: Run the tests and verify the expected failure**

Run: `module load R/4.4.2` followed by `Rscript analysis_v2/tests/testthat.R`

Expected: FAIL because `discover_expression_files()` and `validate_expression_metadata()` do not exist.

- [ ] **Step 3: Add the exact configuration and metadata contract**

Set `analysis.yml` values to:

```yaml
project: TCGA-BRCA
expression_root: GDCdata/TCGA-BRCA/Transcriptome_Profiling/Gene_Expression_Quantification
clinical_patient_file: GDCdata/TCGA-BRCA/Clinical/Clinical_Supplement/8162d394-8b64-4da2-9f5b-d164c54b9608/nationwidechildrens.org_clinical_patient_brca.txt
depmap_root: DepMapCrispr
random_seed: 20260819
expression_filter:
  min_count: 10
  min_fraction: 0.10
differential_expression:
  alpha: 0.05
  lfc_threshold: 1.0
  shrink_type: apeglm
  covariate_min_complete: 0.90
depmap:
  dependency_threshold: -0.5
  dependency_fraction: 0.50
  selectivity_delta: -0.2
```

Set `tractability_evidence.csv` to the header-only contract:

```csv
ensembl_id,gene_symbol,tractability_status,source_name,source_version,source_record_id,evidence_note
```

Append these exact generated-output rules to `.gitignore` while leaving frozen metadata and source manifests trackable:

```gitignore
analysis_v2/results/
analysis_v2/plots/
analysis_v2/provenance/run_log.csv
analysis_v2/provenance/session_info.txt
analysis_v2/provenance/package_versions.csv
```

Implement metadata initialization only when `00_verify_inputs.R --initialize-gdc-metadata` is passed. Query `TCGA-BRCA`, `STAR - Counts`, primary tumor and solid-tissue-normal metadata with `TCGAbiolinks::GDCquery()`/`getResults()`, filter the response to local directory UUIDs, and write:

```text
analysis_v2/provenance/gdc_expression_metadata.tsv
analysis_v2/provenance/gdc_expression_metadata.query.json
analysis_v2/provenance/input_checksums.sha256
```

The frozen table must include `file_id`, `file_name`, `sample_id`, `patient_id`, `sample_type`, retrieval timestamp, and query parameters. A normal invocation must refuse network access and require this table.

The test harness loads `testthat`, sources every file under `analysis_v2/R/` in sorted order, and runs `test_dir("analysis_v2/tests/testthat", reporter = "summary")`. Stage 00 also checks every Tech Stack package and writes its installed version to `analysis_v2/provenance/package_versions.csv`.

- [ ] **Step 4: Run focused and full contract tests**

Run: `Rscript analysis_v2/tests/testthat.R`

Expected: PASS with no network request; tests inject a local `query_fn` fixture.

- [ ] **Step 5: Initialize and verify the real frozen manifest**

Run once: `Rscript analysis_v2/scripts/00_verify_inputs.R --initialize-gdc-metadata`

Then run offline: `Rscript analysis_v2/scripts/00_verify_inputs.R`

Expected: 1,224 unique local files map one-to-one to 1,224 metadata rows; every required input receives a SHA-256 checksum; the offline run reports no API call.

- [ ] **Step 6: Commit the input contract**

```bash
git add .gitignore analysis_v2/config analysis_v2/R/io_contracts.R analysis_v2/scripts/00_verify_inputs.R analysis_v2/provenance analysis_v2/tests
git commit -m "feat: freeze corrected analysis inputs"
```

### Task 2: Offline STAR-count ingestion with stable gene identifiers

**Files:**
- Create: `analysis_v2/R/expression_io.R`
- Create: `analysis_v2/tests/testthat/test-expression-io.R`

**Interfaces:**
- Consumes: frozen metadata from Task 1.
- Produces: `read_star_count_file(path) -> data.frame(ensembl_id, ensembl_id_versioned, gene_symbol, gene_type, count)`
- Produces: `assemble_count_matrix(metadata) -> list(counts, gene_annotation)`
- Produces: `strip_ensembl_version(x) -> character`

- [ ] **Step 1: Write a failing STAR-count fixture test**

```r
test_that("STAR parser removes summary rows and preserves duplicate symbols", {
  path <- tempfile(fileext = ".tsv")
  writeLines(c(
    "# gene-model: GENCODE v36",
    "gene_id\tgene_name\tgene_type\tunstranded",
    "N_unmapped\t\t\t10",
    "ENSG000001.2\tDUP\tprotein_coding\t11",
    "ENSG000002.7\tDUP\tlncRNA\t13"
  ), path)
  x <- read_star_count_file(path)
  expect_identical(x$ensembl_id, c("ENSG000001", "ENSG000002"))
  expect_identical(x$gene_symbol, c("DUP", "DUP"))
  expect_identical(x$count, c(11L, 13L))
})
```

- [ ] **Step 2: Run the fixture test and verify failure**

Run: `Rscript analysis_v2/tests/testthat.R`

Expected: FAIL because `read_star_count_file()` is undefined.

- [ ] **Step 3: Implement strict parsing and matrix assembly**

Use `data.table::fread(skip = 1)`; retain rows matching `^ENSG`; require integer, nonnegative `unstranded` values; and fail if version-stripped Ensembl IDs duplicate within a file. During assembly, require identical ordered Ensembl IDs and annotations in every file, and set matrix column names to `sample_id` from the frozen manifest.

The stage writes `analysis_v2/results/intermediate/counts.rds` and `gene_annotation.csv`, each accompanied by a checksum. It never invokes `GDCprepare()` or `GDCdownload()`.

- [ ] **Step 4: Test malformed counts and inconsistent gene order**

Add cases for a negative count, duplicated Ensembl ID, missing `unstranded`, and two files with different gene orders. Each must fail with the file path and violated invariant in the message.

Run: `Rscript analysis_v2/tests/testthat.R`

Expected: PASS.

- [ ] **Step 5: Commit offline expression ingestion**

```bash
git add analysis_v2/R/expression_io.R analysis_v2/tests/testthat/test-expression-io.R
git commit -m "feat: load frozen STAR counts offline"
```

### Task 3: Clinical normalization and authoritative cohort manifest

**Files:**
- Create: `analysis_v2/R/cohort.R`
- Create: `analysis_v2/scripts/01_build_cohort.R`
- Create: `analysis_v2/tests/testthat/test-cohort.R`

**Interfaces:**
- Consumes: expression metadata and count columns from Tasks 1-2.
- Produces: `read_bcr_patient(path) -> data.frame`
- Produces: `normalize_receptor_status(x) -> factor` with levels `Negative`, `Positive`, `Ambiguous`, `Missing`
- Produces: `classify_tumor(er, pr, her2) -> character`
- Produces: `build_sample_manifest(metadata, clinical) -> data.frame`
- Produces: `select_one_per_patient_group(manifest) -> data.frame`

- [ ] **Step 1: Write the normal-independence regression test**

```r
test_that("normal tissue stays normal regardless of tumor receptor fields", {
  metadata <- data.frame(
    file_id = c("n1", "n2", "t1"),
    sample_id = c("TCGA-AA-0001-11A", "TCGA-AA-0002-11A", "TCGA-AA-0003-01A"),
    patient_id = c("TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003"),
    sample_type = c("Solid Tissue Normal", "Solid Tissue Normal", "Primary Tumor")
  )
  clinical <- data.frame(
    patient_id = metadata$patient_id,
    er_status = c("Positive", NA, "Negative"),
    pr_status = c("Positive", NA, "Negative"),
    her2_status = c("Positive", NA, "Negative")
  )
  manifest <- build_sample_manifest(metadata, clinical)
  expect_identical(manifest$analysis_group, c("Normal", "Normal", "TNBC"))
})
```

- [ ] **Step 2: Verify the regression test fails**

Run: `Rscript analysis_v2/tests/testthat.R`

Expected: FAIL because cohort functions are undefined.

- [ ] **Step 3: Implement receptor and sample classification**

Parse only clinical rows whose patient barcode starts with `TCGA-`, discarding the duplicated descriptive header and `CDE_ID:` row. Classify a primary tumor as `TNBC` only when ER, PR, and HER2 are all `Negative`; classify a receptor-defined other tumor when at least one receptor is `Positive`; otherwise mark it `Unknown_tumor`. Assign `Normal` from sample type before any clinical join and never mutate that assignment from receptor values.

For duplicated samples within patient and group, select the lexicographically smallest full sample barcode and record `deduplication_rank` and `exclusion_reason` for every excluded row.

- [ ] **Step 4: Test ambiguity, pairing, and deterministic deduplication**

Add assertions that equivocal/missing tumors do not enter TNBC, normals remain eligible with missing clinical rows, and shuffled duplicate inputs select the same barcode. Confirm pairing is calculated after one-per-group selection.

Run: `Rscript analysis_v2/tests/testthat.R`

Expected: PASS.

- [ ] **Step 5: Build and audit the real cohort**

Run: `Rscript analysis_v2/scripts/01_build_cohort.R`

Expected outputs:

```text
analysis_v2/results/sample_manifest_all.csv
analysis_v2/results/sample_manifest_included.csv
analysis_v2/results/sample_flow.csv
analysis_v2/results/intermediate/counts.rds
analysis_v2/results/gene_annotation.csv
```

Check that all solid-tissue normals are either included or have only an explicit duplicate/input-integrity exclusion, and that no normal has a receptor-based exclusion.

- [ ] **Step 6: Commit cohort construction**

```bash
git add analysis_v2/R/cohort.R analysis_v2/scripts/01_build_cohort.R analysis_v2/tests/testthat/test-cohort.R
git commit -m "fix: rebuild TNBC and normal cohorts"
```

### Task 4: Primary and sensitivity differential-expression models

**Files:**
- Create: `analysis_v2/R/differential_expression.R`
- Create: `analysis_v2/scripts/02_run_differential_expression.R`
- Create: `analysis_v2/tests/testthat/test-differential-expression.R`

**Interfaces:**
- Consumes: counts, gene annotation, and included sample manifest from Task 3.
- Produces: `filter_counts_for_contrast(counts, sample_ids, min_count, min_fraction) -> matrix`
- Produces: `assess_covariates(col_data, candidates, min_complete) -> data.frame`
- Produces: `run_deseq_contrast(counts, col_data, numerator, denominator, covariates, paired, shrink_type) -> data.frame`
- Produces: `join_gene_annotation(de, genes) -> data.frame`

- [ ] **Step 1: Write failing filter, design, and identifier tests**

```r
test_that("filtering uses only samples in the requested contrast", {
  counts <- matrix(c(10, 10, 0, 0, 999, 999), nrow = 1,
                   dimnames = list("ENSG1", c("a", "b", "outside", "outside2", "x", "y")))
  kept <- filter_counts_for_contrast(counts, c("a", "b"), 10, 0.5)
  expect_identical(rownames(kept), "ENSG1")
  expect_identical(colnames(kept), c("a", "b"))
})

test_that("gene joins retain Ensembl primary keys and duplicate symbols", {
  de <- data.frame(ensembl_id = c("ENSG1", "ENSG2"), padj = c(0.01, 0.02))
  genes <- data.frame(ensembl_id = c("ENSG1", "ENSG2"), gene_symbol = c("DUP", "DUP"))
  joined <- join_gene_annotation(de, genes)
  expect_identical(joined$ensembl_id, c("ENSG1", "ENSG2"))
  expect_identical(joined$gene_symbol, c("DUP", "DUP"))
})
```

- [ ] **Step 2: Run tests and verify failure**

Run: `Rscript analysis_v2/tests/testthat.R`

Expected: FAIL because DE functions are undefined.

- [ ] **Step 3: Implement model gates and DESeq2 output contract**

Require at least two samples per group, a full-rank design matrix, and no missing design values. A candidate covariate is eligible only at or above 90% completeness, with variation in both groups, no level confined to one group, and a full-rank combined design. Write every accepted/rejected decision and reason to `covariate_eligibility.csv`.

Return one row per Ensembl ID with:

```text
ensembl_id, gene_symbol, gene_type, base_mean,
log2fc_raw, lfc_se_raw, stat, pvalue, padj,
log2fc_shrunk, lfc_se_shrunk, contrast, n_numerator, n_denominator
```

Use `DESeq2::lfcShrink(..., type = "apeglm")` on the named numerator-versus-denominator coefficient and seed stochastic operations with `20260819`.

- [ ] **Step 4: Implement the three explicit contrasts**

Write:

```text
analysis_v2/results/de_primary_tnbc_vs_normal.csv
analysis_v2/results/de_paired_tnbc_vs_normal.csv
analysis_v2/results/de_specificity_tnbc_vs_other.csv
analysis_v2/results/covariate_eligibility.csv
analysis_v2/results/de_contrast_summary.csv
```

The paired model uses only patients with one included TNBC and one included normal and a design of `~ patient_id + group`. The specificity model excludes `Unknown_tumor` samples.

- [ ] **Step 5: Run unit tests and real models**

Run: `Rscript analysis_v2/tests/testthat.R`

Then: `Rscript analysis_v2/scripts/02_run_differential_expression.R`

Expected: all three result tables have unique Ensembl IDs and complete contrast metadata; the paired model reports the exact number of pairs; no `.1`-style symbols appear as generated identifiers.

- [ ] **Step 6: Commit differential expression**

```bash
git add analysis_v2/R/differential_expression.R analysis_v2/scripts/02_run_differential_expression.R analysis_v2/tests/testthat/test-differential-expression.R
git commit -m "feat: add corrected TNBC expression models"
```

### Task 5: Universe-aware ORA and ranked enrichment

**Files:**
- Create: `analysis_v2/R/enrichment.R`
- Create: `analysis_v2/scripts/03_run_enrichment.R`
- Create: `analysis_v2/tests/testthat/test-enrichment.R`

**Interfaces:**
- Consumes: DE tables from Task 4.
- Produces: `prepare_ora_inputs(de, term2gene, direction, alpha, lfc_threshold) -> list(genes, universe)`
- Produces: `map_tested_genes(de, org_db) -> list(mapped, audit)`
- Produces: `run_ora(de, term2gene, direction, alpha, lfc_threshold) -> data.frame`
- Produces: `run_ranked_enrichment(de, term2gene, seed) -> data.frame`
- Produces: `build_kegg_legacy_mapping() -> data.frame(term, ensembl_id)`

- [ ] **Step 1: Write a failing universe regression test**

```r
test_that("ORA background contains only tested and mapped genes", {
  de <- data.frame(
    ensembl_id = c("E1", "E2", "E3"),
    padj = c(0.01, 0.20, NA),
    log2fc_shrunk = c(2, 0, 1)
  )
  terms <- data.frame(term = c("P", "P", "P"), ensembl_id = c("E1", "E2", "E9"))
  inputs <- prepare_ora_inputs(de, terms, "up", 0.05, 1)
  expect_setequal(inputs$universe, c("E1", "E2"))
  expect_identical(inputs$genes, "E1")
})
```

- [ ] **Step 2: Run the test and verify failure**

Run: `Rscript analysis_v2/tests/testthat.R`

Expected: FAIL because `prepare_ora_inputs()` is undefined.

- [ ] **Step 3: Implement local, version-recorded mappings**

Use `org.Hs.eg.db`/`GO.db` for GO Biological Process and the installed `msigdbr` `C2:CP:KEGG_LEGACY` snapshot for KEGG-compatible pathways. Do not call KEGGREST or any online enrichment endpoint. Record package/database versions and mapping counts.

For ORA, pass the tested-and-mapped universe explicitly. For ranked enrichment, use unique Ensembl IDs ranked by the signed DE statistic, set the seed, and report duplicate/missing mapping losses.

- [ ] **Step 4: Write directional and empty-set tests**

Test up and down selection separately. An empty significant set must return a zero-row table plus an audit row with `status = "no_significant_genes"`; it must not error or create a misleading plot.

Run: `Rscript analysis_v2/tests/testthat.R`

Expected: PASS.

- [ ] **Step 5: Run enrichment and verify audit outputs**

Run: `Rscript analysis_v2/scripts/03_run_enrichment.R`

Expected outputs include GO and KEGG-legacy ORA for both directions, ranked enrichment for the primary contrast, `mapping_audit.csv`, and `enrichment_run_metadata.csv`. Every ORA row must report foreground and universe sizes.

- [ ] **Step 6: Commit enrichment**

```bash
git add analysis_v2/R/enrichment.R analysis_v2/scripts/03_run_enrichment.R analysis_v2/tests/testthat/test-enrichment.R
git commit -m "feat: add universe-aware pathway analysis"
```

### Task 6: Multi-model DepMap integration

**Files:**
- Create: `analysis_v2/R/depmap.R`
- Create: `analysis_v2/scripts/04_integrate_depmap.R`
- Create: `analysis_v2/tests/testthat/test-depmap.R`

**Interfaces:**
- Consumes: primary DE results and local DepMap CSVs.
- Produces: `classify_depmap_models(models) -> data.frame(model_id, model_name, model_group, source_label)`
- Produces: `extract_depmap_effects(effect_matrix, model_manifest) -> data.frame`
- Produces: `summarize_depmap_by_gene(long_effects, threshold) -> data.frame`
- Produces: `annotate_common_essentials(summary, common_table) -> data.frame`

- [ ] **Step 1: Write failing model-group and summary tests**

```r
test_that("only explicit breast subtype labels enter comparator groups", {
  models <- data.frame(
    ModelID = c("T", "E", "U", "O"),
    CellLineName = c("TN", "ER", "Unknown", "Ovary"),
    OncotreeLineage = c("Breast", "Breast", "Breast", "Ovary/Fallopian Tube"),
    ModelSubtypeFeatures = c("basal_A TNBC", "luminal ER+", NA, "TNBC")
  )
  groups <- classify_depmap_models(models)
  expect_identical(groups$model_id, c("T", "E"))
  expect_identical(groups$model_group, c("TNBC", "Non_TNBC_breast"))
})

test_that("dependency summaries retain panel size and selectivity direction", {
  x <- data.frame(gene_symbol = rep("G", 4), model_group = rep(c("TNBC", "Non_TNBC_breast"), each = 2), effect = c(-1, -0.8, -0.2, 0))
  out <- summarize_depmap_by_gene(x, -0.5)
  expect_equal(out$tnbc_median_effect, -0.9)
  expect_equal(out$non_tnbc_median_effect, -0.1)
  expect_equal(out$tnbc_minus_non_tnbc, -0.8)
})
```

- [ ] **Step 2: Run tests and verify failure**

Run: `Rscript analysis_v2/tests/testthat.R`

Expected: FAIL because the new DepMap functions are undefined.

- [ ] **Step 3: Implement transparent model grouping**

Use `ModelSubtypeFeatures` literally: breast models containing `TNBC` enter the TNBC panel; breast models not containing `TNBC` but explicitly containing `ER+` or `HER2+` enter the non-TNBC comparator. Exclude unlabeled breast models and export their reasons. Preserve the complete source label in the manifest. Confirm MDA-MB-231 resolves uniquely to `ACH-000768`.

- [ ] **Step 4: Implement gene-effect and common-essential summaries**

Normalize DepMap headers by removing the trailing numeric Entrez annotation only. Produce per-gene TNBC and non-TNBC model counts, median effects, interquartile ranges, fractions at or below `-0.5`, TNBC-minus-non-TNBC difference, MDA-MB-231 effect, and common-essential status. Do not treat `ScreenGeneDependency.csv` as model-level probabilities.

- [ ] **Step 5: Run tests and the real integration**

Run: `Rscript analysis_v2/tests/testthat.R`

Then: `Rscript analysis_v2/scripts/04_integrate_depmap.R`

Expected outputs:

```text
analysis_v2/results/depmap_model_manifest.csv
analysis_v2/results/depmap_model_exclusions.csv
analysis_v2/results/depmap_gene_dependency_summary.csv
analysis_v2/results/de_primary_with_depmap.csv
analysis_v2/results/depmap_run_metadata.csv
```

The report must label the DepMap release unresolved unless an exact release identifier exists in a local source file.

- [ ] **Step 6: Commit DepMap integration**

```bash
git add analysis_v2/R/depmap.R analysis_v2/scripts/04_integrate_depmap.R analysis_v2/tests/testthat/test-depmap.R
git commit -m "feat: integrate multi-model DepMap evidence"
```

### Task 7: Evidence classes, figures, and scientific report

**Files:**
- Create: `analysis_v2/R/evidence.R`
- Create: `analysis_v2/R/reporting.R`
- Create: `analysis_v2/scripts/05_build_report.R`
- Create: `analysis_v2/corrected_analysis_report.Rmd`
- Create: `analysis_v2/tests/testthat/test-evidence-reporting.R`

**Interfaces:**
- Consumes: primary, paired, specificity, enrichment, and DepMap tables.
- Produces: `assemble_candidate_evidence(primary, paired, specificity, depmap, tractability, config) -> data.frame`
- Produces: `assign_evidence_class(x, config) -> character`
- Produces: `validate_final_outputs(root) -> data.frame`
- Produces: `write_output_manifest(paths, output_file) -> data.frame`
- Produces: `render_corrected_report(output_file) -> character`

- [ ] **Step 1: Write a failing no-composite-score test**

```r
test_that("candidate evidence is classified without an arithmetic score", {
  x <- data.frame(
    primary_hit = c(TRUE, TRUE, TRUE),
    paired_concordant = c(TRUE, FALSE, FALSE),
    specificity_concordant = c(TRUE, TRUE, FALSE),
    tnbc_dependency = c(TRUE, TRUE, FALSE),
    depmap_selective = c(TRUE, FALSE, FALSE),
    common_essential = c(FALSE, FALSE, FALSE)
  )
  out <- assign_evidence_class(x, list())
  expect_identical(out, c("A_multistream", "B_expression_dependency", "C_expression_only"))
  expect_false(any(grepl("score", names(x), ignore.case = TRUE)))
})
```

- [ ] **Step 2: Run tests and verify failure**

Run: `Rscript analysis_v2/tests/testthat.R`

Expected: FAIL because evidence functions are undefined.

- [ ] **Step 3: Implement explicit evidence flags and classes**

Define `primary_hit` as primary FDR below 0.05 and shrunken log2FC above 1. Define paired and specificity concordance by positive effect direction, recording unavailable results separately. Define `tnbc_dependency` as TNBC median effect at or below `-0.5` and at least half of measured TNBC models at or below `-0.5`; define DepMap selectivity as TNBC-minus-non-TNBC at or below `-0.2`.

Assign:

- `A_multistream`: primary hit, both sensitivity directions concordant, TNBC dependency, DepMap selectivity, and not a common essential.
- `B_expression_dependency`: primary hit, TNBC dependency, and at least one concordant sensitivity direction.
- `C_expression_only`: primary hit not meeting A or B.
- `not_candidate`: not a primary hit.

Keep every underlying field. `tractability_evidence.csv` contains only its header until a versioned source is supplied, so the report must say `not assessed`, not `not druggable`.

- [ ] **Step 4: Implement report checks and figures**

Generate PCA, primary volcano, paired-versus-primary concordance, specificity-versus-primary concordance, enrichment dot plots when nonempty, and TNBC-versus-non-TNBC dependency plots. Give every plot sample sizes, thresholds, and contrast direction in its caption.

The R Markdown report must include cohort flow, covariate decisions, primary and sensitivity results, mapping losses, enrichment universes, DepMap panel composition, candidate evidence classes, input provenance, warnings, and limitations from the spec. Stage 05 also writes `analysis_v2/provenance/output_manifest.sha256`, `warnings.csv`, and `session_info.txt`; validation fails if a file claimed by the report is absent from the output manifest.

- [ ] **Step 5: Run report tests and render**

Run: `Rscript analysis_v2/tests/testthat.R`

Then: `Rscript analysis_v2/scripts/05_build_report.R`

Expected: `analysis_v2/results/candidate_evidence.csv`, nonzero PNG files under `analysis_v2/plots/`, `analysis_v2/results/corrected_analysis_report.html`, and the three run-provenance records. Open each representative PNG and the HTML report for visual inspection; reject clipped labels, unreadable legends, missing captions, or claims that exceed target nomination.

- [ ] **Step 6: Commit reporting**

```bash
git add analysis_v2/R/evidence.R analysis_v2/R/reporting.R analysis_v2/scripts/05_build_report.R analysis_v2/corrected_analysis_report.Rmd analysis_v2/tests/testthat/test-evidence-reporting.R
git commit -m "feat: report separated TNBC target evidence"
```

### Task 8: End-to-end runner, documentation, and final verification

**Files:**
- Create: `analysis_v2/run_all.R`
- Create: `analysis_v2/README.md`
- Create: `analysis_v2/provenance/source_checksums.sha256`
- Modify: `AGENTS.md`
- Modify: `PROJECT_RESTART_STATUS.md`

**Interfaces:**
- Consumes: all stage scripts and contracts from Tasks 1-7.
- Produces: a fail-fast, offline rerun and a documented operator handoff.

- [ ] **Step 1: Write the runner with stage boundaries**

`run_all.R` must call stages 00-05 in order with `system2("Rscript", stage)`, stop on the first nonzero status, and write start/end time, status, and command for each stage to `analysis_v2/provenance/run_log.csv`. It must not pass the metadata-initialization flag.

- [ ] **Step 2: Document exact operation and recovery commands**

`analysis_v2/README.md` must document:

```bash
module load R/4.4.2
Rscript analysis_v2/tests/testthat.R
Rscript analysis_v2/scripts/00_verify_inputs.R
Rscript analysis_v2/run_all.R
```

Explain how to resume from a numbered stage, which outputs are generated, why the one-time metadata file is required, and why legacy outputs are not comparable without checking cohort and model definitions.

- [ ] **Step 3: Run static and unit verification**

Run:

```bash
Rscript -e 'files <- c(Sys.glob("analysis_v2/R/*.R"), Sys.glob("analysis_v2/scripts/*.R")); invisible(lapply(files, parse))'
Rscript analysis_v2/tests/testthat.R
```

Expected: every script parses and all tests pass.

- [ ] **Step 4: Run the complete frozen-data workflow**

Run: `Rscript analysis_v2/run_all.R`

Expected: all six stages report success, required result tables are nonempty where scientifically expected, checksums match, and the report records the exact sample and gene counts. If a scientifically valid result set is empty, the stage succeeds only with an explicit audit status and no misleading plot.

- [ ] **Step 5: Verify legacy preservation and provenance**

Run:

```bash
awk '$2 != "DepMapCrispr/DepMapCrispr.R"' provenance/baseline_2025-09-22.sha256 | sha256sum -c -
sha256sum -c provenance/restart_2026-08-19.sha256
sha256sum -c analysis_v2/provenance/source_checksums.sha256
git status --short -- .
```

Expected: preserved baseline and restart entries pass. In `05_build_report.R`, construct the sorted source path vector with `list.files(..., recursive = TRUE, full.names = TRUE)` over `analysis_v2/R/`, `analysis_v2/scripts/`, `analysis_v2/config/`, and `analysis_v2/tests/`, append `analysis_v2/run_all.R`, the report template, the approved spec, and this plan, then call Task 1's `write_checksum_manifest()` to generate `source_checksums.sha256`. Exclude generated results, plots, run logs, and package/session records from the source path vector.

- [ ] **Step 6: Update contributor and restart documentation**

Add the corrected stage commands and directory map to `AGENTS.md`. Update `PROJECT_RESTART_STATUS.md` with the corrected cohort counts, test count, completed stages, output locations, unresolved release labels, and remaining scientific limitations.

- [ ] **Step 7: Commit the completed corrected workflow**

```bash
git add AGENTS.md PROJECT_RESTART_STATUS.md analysis_v2 docs/superpowers/specs/2026-08-19-corrected-tnbc-analysis-design.md docs/superpowers/plans/2026-08-19-corrected-tnbc-analysis.md
git commit -m "feat: add reproducible TNBC correction workflow"
```
