#!/usr/bin/env Rscript

required_packages <- c("data.table", "dplyr", "janitor", "ggplot2", "ggrepel")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0) {
  stop("Install required packages before running: ", paste(missing_packages, collapse = ", "))
}

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_argument) != 1) {
  stop("Run this integration with Rscript")
}
script_file <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
script_dir <- dirname(script_file)

source(file.path(script_dir, "R", "depmap_integration.R"))

invisible(run_depmap_integration(
  input_dir = script_dir,
  output_dir = file.path(script_dir, "depmap_integration_outputs")
))
