#!/usr/bin/env Rscript

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_argument) != 1) {
  stop("Run this stage with Rscript", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
analysis_root <- file.path(project_root, "analysis_v2")

source(file.path(analysis_root, "R", "io_contracts.R"))

arguments <- commandArgs(trailingOnly = TRUE)
allowed_arguments <- "--initialize-gdc-metadata"
unexpected_arguments <- setdiff(arguments, allowed_arguments)
if (length(unexpected_arguments) > 0) {
  stop("Unknown arguments: ", paste(unexpected_arguments, collapse = ", "), call. = FALSE)
}
initialize_metadata <- allowed_arguments %in% arguments

required_packages <- c(
  "data.table", "dplyr", "yaml", "digest", "jsonlite", "DESeq2", "apeglm",
  "clusterProfiler", "AnnotationDbi", "org.Hs.eg.db", "GO.db", "msigdbr",
  "fgsea", "ggplot2", "ggrepel", "rmarkdown", "testthat", "TCGAbiolinks"
)
available <- vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
if (any(!available)) {
  stop(
    "Missing required R packages: ",
    paste(required_packages[!available], collapse = ", "),
    call. = FALSE
  )
}

config <- read_analysis_config(file.path(analysis_root, "config", "analysis.yml"))
expression_root <- file.path(project_root, config$expression_root)
local_files <- discover_expression_files(expression_root)
metadata_path <- file.path(analysis_root, "provenance", "gdc_expression_metadata.tsv")

metadata <- load_or_initialize_gdc_metadata(
  local_files = local_files,
  metadata_path = metadata_path,
  initialize = initialize_metadata,
  query_fn = function() query_gdc_expression_metadata(config)
)

if (initialize_metadata) {
  query_record <- list(
    lookup_kind = "metadata_only",
    project = config$project,
    data_category = "Transcriptome Profiling",
    data_type = "Gene Expression Quantification",
    workflow_type = "STAR - Counts",
    sample_type = c("Primary Tumor", "Solid Tissue Normal"),
    local_file_count = nrow(local_files),
    downloaded_data = FALSE,
    retrieved_at_utc = unique(metadata$retrieved_at_utc)
  )
  jsonlite::write_json(
    query_record,
    file.path(analysis_root, "provenance", "gdc_expression_metadata.query.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )
}

depmap_files <- file.path(
  project_root,
  config$depmap_root,
  c(
    "CRISPRGeneEffect.csv", "CRISPRInferredCommonEssentials.csv", "Gene.csv",
    "Model.csv", "ScreenGeneDependency.csv", "ScreenGeneEffect.csv"
  )
)
required_inputs <- c(
  local_files$path,
  file.path(project_root, config$clinical_patient_file),
  depmap_files,
  metadata_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing required inputs: ", paste(missing_inputs, collapse = ", "), call. = FALSE)
}

package_versions <- data.frame(
  package = required_packages,
  version = vapply(
    required_packages,
    function(package) as.character(utils::packageVersion(package)),
    FUN.VALUE = character(1)
  ),
  stringsAsFactors = FALSE
)
data.table::fwrite(
  package_versions,
  file.path(analysis_root, "provenance", "package_versions.csv")
)
invisible(write_checksum_manifest(
  required_inputs,
  file.path(analysis_root, "provenance", "input_checksums.sha256")
))

message("Verified ", nrow(local_files), " frozen expression files and required clinical/DepMap inputs.")
message("GDC sample metadata rows: ", nrow(metadata))
