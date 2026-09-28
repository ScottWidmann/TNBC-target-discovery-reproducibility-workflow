#!/usr/bin/env Rscript

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_argument) != 1) {
  stop("Run this stage with Rscript", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
analysis_root <- file.path(project_root, "analysis_v2")

source(file.path(analysis_root, "R", "io_contracts.R"))
source(file.path(analysis_root, "R", "expression_io.R"))
source(file.path(analysis_root, "R", "cohort.R"))

config <- read_analysis_config(file.path(analysis_root, "config", "analysis.yml"))
local_files <- discover_expression_files(file.path(project_root, config$expression_root))
metadata <- load_or_initialize_gdc_metadata(
  local_files,
  file.path(analysis_root, "provenance", "gdc_expression_metadata.tsv"),
  initialize = FALSE,
  query_fn = function() stop("Stage 01 must not query GDC", call. = FALSE)
)
metadata$path <- local_files$path[match(metadata$file_id, local_files$file_id)]
if (!identical(metadata$file_name, local_files$file_name[match(metadata$file_id, local_files$file_id)])) {
  stop("Frozen GDC metadata file names do not match local files", call. = FALSE)
}

clinical <- read_bcr_patient(file.path(project_root, config$clinical_patient_file))
manifest <- build_sample_manifest(metadata, clinical)
manifest <- select_one_per_patient_group(manifest)
included <- manifest[manifest$selected_for_analysis, , drop = FALSE]
if (nrow(included) == 0) {
  stop("Cohort construction selected no samples", call. = FALSE)
}
if (!all(c("TNBC", "Normal", "Other_receptor_defined") %in% included$analysis_group)) {
  stop("Cohort is missing a required analysis group", call. = FALSE)
}

results_dir <- file.path(analysis_root, "results")
intermediate_dir <- file.path(results_dir, "intermediate")
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

data.table::fwrite(manifest, file.path(results_dir, "sample_manifest_all.csv"))
data.table::fwrite(included, file.path(results_dir, "sample_manifest_included.csv"))

sample_flow <- data.table::as.data.table(manifest)[, .(
  total_files = .N,
  eligible_files = sum(eligible),
  selected_samples = sum(selected_for_analysis),
  excluded_files = sum(!selected_for_analysis)
), by = .(sample_type, analysis_group, exclusion_reason)]
data.table::setorder(sample_flow, sample_type, analysis_group, exclusion_reason)
data.table::fwrite(sample_flow, file.path(results_dir, "sample_flow.csv"))

assembled <- assemble_count_matrix(included[, c("file_id", "path", "sample_id")])
if (!identical(colnames(assembled$counts), included$sample_id)) {
  stop("Count matrix columns do not match included sample manifest", call. = FALSE)
}
saveRDS(
  assembled$counts,
  file.path(intermediate_dir, "counts.rds"),
  compress = FALSE
)
data.table::fwrite(
  assembled$gene_annotation,
  file.path(results_dir, "gene_annotation.csv")
)

cohort_outputs <- c(
  file.path(results_dir, "sample_manifest_all.csv"),
  file.path(results_dir, "sample_manifest_included.csv"),
  file.path(results_dir, "sample_flow.csv"),
  file.path(results_dir, "gene_annotation.csv"),
  file.path(intermediate_dir, "counts.rds")
)
invisible(write_checksum_manifest(
  cohort_outputs,
  file.path(analysis_root, "provenance", "cohort_outputs.sha256")
))

message("Selected samples by analysis group:")
print(table(included$analysis_group))
message("Paired TNBC-normal patients: ", length(unique(
  included$patient_id[included$paired_tnbc_normal]
)))
message("Count matrix: ", nrow(assembled$counts), " genes x ", ncol(assembled$counts), " samples")
