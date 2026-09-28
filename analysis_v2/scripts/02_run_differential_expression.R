#!/usr/bin/env Rscript

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_argument) != 1) {
  stop("Run this stage with Rscript", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
analysis_root <- file.path(project_root, "analysis_v2")

source(file.path(analysis_root, "R", "io_contracts.R"))
source(file.path(analysis_root, "R", "differential_expression.R"))

config <- read_analysis_config(file.path(analysis_root, "config", "analysis.yml"))
results_dir <- file.path(analysis_root, "results")
counts <- readRDS(file.path(results_dir, "intermediate", "counts.rds"))
manifest <- data.table::fread(
  file.path(results_dir, "sample_manifest_included.csv"),
  data.table = FALSE
)
genes <- data.table::fread(
  file.path(results_dir, "gene_annotation.csv"),
  data.table = FALSE
)
if (!identical(colnames(counts), manifest$sample_id)) {
  stop("Count matrix columns do not match the included sample manifest", call. = FALSE)
}

set.seed(config$random_seed)
filter_config <- config$expression_filter
de_config <- config$differential_expression

run_declared_contrast <- function(
  name,
  numerator,
  denominator,
  paired = FALSE
) {
  selected <- manifest$analysis_group %in% c(numerator, denominator)
  if (paired) selected <- selected & manifest$paired_tnbc_normal
  col_data <- manifest[selected, , drop = FALSE]
  col_data$group <- col_data$analysis_group
  covariates <- character()
  if (!paired) {
    audit <- assess_covariates(
      col_data,
      de_config$candidate_covariates,
      de_config$covariate_min_complete
    )
    audit$contrast <- name
    covariates <- audit$covariate[audit$eligible]
  } else {
    audit <- data.frame(
      covariate = "patient_id",
      completeness = 1,
      eligible = TRUE,
      reason = "paired_design_block",
      contrast = name,
      stringsAsFactors = FALSE
    )
  }
  filtered_counts <- filter_counts_for_contrast(
    counts,
    col_data$sample_id,
    filter_config$min_count,
    filter_config$min_fraction
  )
  result <- run_deseq_contrast(
    filtered_counts,
    col_data,
    numerator = numerator,
    denominator = denominator,
    covariates = covariates,
    paired = paired,
    shrink_type = de_config$shrink_type
  )
  result <- join_gene_annotation(result, genes)
  summary <- data.frame(
    contrast = name,
    numerator = numerator,
    denominator = denominator,
    paired = paired,
    covariates = if (length(covariates) == 0) "none" else paste(covariates, collapse = ";"),
    numerator_samples = unique(result$n_numerator),
    denominator_samples = unique(result$n_denominator),
    tested_genes = nrow(result),
    fdr_below_alpha = sum(result$padj < de_config$alpha, na.rm = TRUE),
    upregulated = sum(
      result$padj < de_config$alpha &
        result$log2fc_shrunk > de_config$lfc_threshold,
      na.rm = TRUE
    ),
    downregulated = sum(
      result$padj < de_config$alpha &
        result$log2fc_shrunk < -de_config$lfc_threshold,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
  list(result = result, audit = audit, summary = summary)
}

primary <- run_declared_contrast(
  "primary_tnbc_vs_normal",
  numerator = "TNBC",
  denominator = "Normal"
)
paired <- run_declared_contrast(
  "paired_tnbc_vs_normal",
  numerator = "TNBC",
  denominator = "Normal",
  paired = TRUE
)
specificity <- run_declared_contrast(
  "specificity_tnbc_vs_other",
  numerator = "TNBC",
  denominator = "Other_receptor_defined"
)

output_paths <- c(
  primary = file.path(results_dir, "de_primary_tnbc_vs_normal.csv"),
  paired = file.path(results_dir, "de_paired_tnbc_vs_normal.csv"),
  specificity = file.path(results_dir, "de_specificity_tnbc_vs_other.csv"),
  covariates = file.path(results_dir, "covariate_eligibility.csv"),
  summary = file.path(results_dir, "de_contrast_summary.csv")
)
data.table::fwrite(primary$result, output_paths[["primary"]])
data.table::fwrite(paired$result, output_paths[["paired"]])
data.table::fwrite(specificity$result, output_paths[["specificity"]])
data.table::fwrite(
  rbind(primary$audit, paired$audit, specificity$audit),
  output_paths[["covariates"]]
)
contrast_summary <- rbind(primary$summary, paired$summary, specificity$summary)
data.table::fwrite(contrast_summary, output_paths[["summary"]])
invisible(write_checksum_manifest(
  unname(output_paths),
  file.path(analysis_root, "provenance", "de_outputs.sha256")
))

print(contrast_summary)
