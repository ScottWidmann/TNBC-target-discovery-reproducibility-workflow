#!/usr/bin/env Rscript

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_argument) != 1) {
  stop("Run this workflow with Rscript", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
analysis_root <- normalizePath(dirname(script_file), mustWork = TRUE)
project_root <- normalizePath(file.path(analysis_root, ".."), mustWork = TRUE)

source(file.path(analysis_root, "R", "runner.R"))
setwd(project_root)
Sys.setenv(TZ = "UTC")

stages <- file.path(
  analysis_root,
  "scripts",
  sprintf("%02d_%s.R", 0:5, c(
    "verify_inputs",
    "build_cohort",
    "run_differential_expression",
    "run_enrichment",
    "integrate_depmap",
    "build_report"
  ))
)
invisible(run_analysis_stages(
  stages,
  file.path(analysis_root, "provenance", "run_log.csv")
))
message("Corrected TNBC workflow completed successfully")
