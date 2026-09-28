# Corrected TNBC Target-Nomination Workflow

`analysis_v2/` is the reproducible correction of the preserved exploratory
workflow. It reads the existing TCGA-BRCA STAR counts, clinical file, and DepMap
CSVs without downloading or replacing them. It never overwrites legacy
`results/`, `plots/`, or scripts.

## Run and Test

From the repository root:

```bash
module load R/4.4.2
Rscript analysis_v2/tests/testthat.R
Rscript analysis_v2/scripts/00_verify_inputs.R
Rscript analysis_v2/run_all.R
```

Stage 00 requires the frozen
`analysis_v2/provenance/gdc_expression_metadata.tsv`. It was created by the one
approved metadata-only GDC lookup that maps local file UUIDs to sample IDs. A
normal verification or full run does not query GDC and fails closed if the
metadata or frozen inputs do not match.

## Stages and Recovery

The runner executes these restartable stages in order:

1. `00_verify_inputs.R`: validate packages, schemas, metadata, and input checksums.
2. `01_build_cohort.R`: assemble counts and the authoritative sample manifest.
3. `02_run_differential_expression.R`: run primary, paired, and specificity models.
4. `03_run_enrichment.R`: run local GO/KEGG ORA and ranked enrichment.
5. `04_integrate_depmap.R`: summarize TNBC and comparator breast models.
6. `05_build_report.R`: assign evidence classes, create plots, and render HTML.

If a stage fails, inspect `analysis_v2/provenance/run_log.csv`, correct the
cause, and rerun that stage directly, for example:

```bash
Rscript analysis_v2/scripts/03_run_enrichment.R
Rscript analysis_v2/scripts/04_integrate_depmap.R
Rscript analysis_v2/scripts/05_build_report.R
```

Each stage validates its upstream contracts. Do not skip a failed stage or
present downstream files from a partial run as current.

## Outputs and Interpretation

Tables are written to `analysis_v2/results/`, figures to `analysis_v2/plots/`,
and checksums, warnings, versions, and run records to
`analysis_v2/provenance/`. The main deliverables are
`candidate_evidence.csv` and `corrected_analysis_report.html`.

Legacy and corrected results are not directly interchangeable: they use
different normal cohorts, covariate rules, sensitivity contrasts, gene-key
contracts, DepMap panels, enrichment universes, and candidate definitions.
Evidence classes preserve individual measurements and use no arithmetic score.
They support target nomination only; tractability remains not assessed until a
versioned curated source is supplied.
