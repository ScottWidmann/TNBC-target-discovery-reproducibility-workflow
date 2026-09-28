# Repository Guidelines

## Project Structure & Module Organization

This repository contains two distinct R workflows. Treat `tnbc_analysis_claude.r`, legacy `results/` and `plots/`, and `DepMapCrispr/DepMapCrispr.R` as preserved exploratory history. The corrected, reproducible workflow is under `analysis_v2/`: reusable functions are in `R/`, numbered stages in `scripts/`, settings in `config/`, tests in `tests/testthat/`, and audit records in `provenance/`. Generated corrected tables and figures belong only in `analysis_v2/results/` and `analysis_v2/plots/`. Frozen TCGA files remain under `GDCdata/`; local DepMap inputs remain under `DepMapCrispr/`.

## Build, Test, and Development Commands

Load the project R version and run from the repository root:

```bash
module load R/4.4.2
Rscript analysis_v2/tests/testthat.R
Rscript analysis_v2/scripts/00_verify_inputs.R
Rscript analysis_v2/run_all.R
```

The test command runs fixture and regression tests. Stage 00 verifies frozen inputs and performs no network query during normal use. `run_all.R` executes stages 00–05 offline and stops at the first failure. To resume, invoke the failed numbered stage and then later stages explicitly. Parse changed scripts with `Rscript -e 'parse(file="path/to/script.R")'`.

## Coding Style & Naming Conventions

Use two-space indentation, UTF-8, spaces instead of tabs, `snake_case`, and `<-` assignment. Prefer explicit package qualification such as `dplyr::filter`. Keep Ensembl IDs as primary gene keys; symbols are annotations and may duplicate. Put thresholds in `analysis_v2/config/analysis.yml`, not inline constants. Never silently substitute inputs, models, or smaller cohorts.

## Testing Guidelines

Add a failing `testthat` fixture before correcting behavior, then run the full suite. Test schema failures, unique-key invariants, unavailable evidence, and deterministic selection. For real-data stages, also verify checksum manifests, nonempty expected outputs, and figure/report readability. A valid empty scientific result needs an explicit status, not a fabricated plot.

## Commit & Pull Request Guidelines

Use concise Conventional Commit subjects, for example `fix: preserve normal tissue eligibility`. Pull requests should state the scientific question, frozen input snapshot, cohort/model changes, commands run, test results, and limitations; include representative figures when visuals change. Do not commit credentials, `.Rhistory`, `.Rproj.user/`, generated results, or redundant large downloads. Do not initialize, commit, or push this directory into its enclosing workspace repository.
