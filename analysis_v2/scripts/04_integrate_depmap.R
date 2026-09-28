#!/usr/bin/env Rscript

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_argument) != 1) {
  stop("Run this stage with Rscript", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
analysis_root <- file.path(project_root, "analysis_v2")

source(file.path(analysis_root, "R", "io_contracts.R"))
source(file.path(analysis_root, "R", "depmap.R"))

config <- read_analysis_config(file.path(analysis_root, "config", "analysis.yml"))
results_dir <- file.path(analysis_root, "results")
depmap_root <- file.path(project_root, config$depmap_root)

model_path <- file.path(depmap_root, "Model.csv")
effect_path <- file.path(depmap_root, "CRISPRGeneEffect.csv")
common_path <- file.path(depmap_root, "CRISPRInferredCommonEssentials.csv")
models <- data.table::fread(model_path, data.table = FALSE)
classified <- classify_depmap_models(models)

effect_ids <- data.table::fread(effect_path, select = 1, data.table = FALSE)
effect_id_column <- names(effect_ids)[[1]]
measured <- classified$model_id %in% as.character(effect_ids[[effect_id_column]])
model_manifest <- classified[measured, , drop = FALSE]
model_manifest$measured_gene_effect <- TRUE

is_breast <- toupper(trimws(models$OncotreeLineage)) == "BREAST"
breast_models <- data.frame(
  model_id = as.character(models$ModelID[is_breast]),
  model_name = as.character(models$CellLineName[is_breast]),
  source_label = as.character(models$ModelSubtypeFeatures[is_breast]),
  stringsAsFactors = FALSE
)
breast_models$classified <- breast_models$model_id %in% classified$model_id
breast_models$measured <- breast_models$model_id %in% model_manifest$model_id
model_exclusions <- breast_models[!breast_models$measured, , drop = FALSE]
model_exclusions$exclusion_reason <- ifelse(
  !model_exclusions$classified,
  "subtype_label_not_explicit_tnbc_er_or_her2",
  "missing_from_gene_effect_matrix"
)

if (!all(c("TNBC", "Non_TNBC_breast") %in% model_manifest$model_group)) {
  stop("Measured DepMap manifest lacks a required comparison group", call. = FALSE)
}
if (sum(toupper(gsub("[^A-Z0-9]", "", model_manifest$model_name)) == "MDAMB231") != 1) {
  stop("MDA-MB-231 does not resolve uniquely in the measured TNBC panel", call. = FALSE)
}

message(
  "Reading frozen DepMap gene-effect matrix for ", nrow(model_manifest),
  " classified breast models"
)
effects <- data.table::fread(effect_path, data.table = FALSE, check.names = FALSE)
long_effects <- extract_depmap_effects(effects, model_manifest)
rm(effects)
invisible(gc())

dependency <- summarize_depmap_by_gene(
  long_effects,
  dependency_threshold = config$depmap$dependency_threshold
)
common_table <- read_common_essential_table(common_path)
dependency <- annotate_common_essentials(dependency, common_table)

primary_de <- data.table::fread(
  file.path(results_dir, "de_primary_tnbc_vs_normal.csv"),
  data.table = FALSE
)
joined <- join_depmap_to_de(primary_de, dependency)

output_paths <- c(
  manifest = file.path(results_dir, "depmap_model_manifest.csv"),
  exclusions = file.path(results_dir, "depmap_model_exclusions.csv"),
  dependency = file.path(results_dir, "depmap_gene_dependency_summary.csv"),
  joined = file.path(results_dir, "de_primary_with_depmap.csv"),
  metadata = file.path(results_dir, "depmap_run_metadata.csv")
)
data.table::fwrite(model_manifest, output_paths[["manifest"]])
data.table::fwrite(model_exclusions, output_paths[["exclusions"]])
data.table::fwrite(dependency, output_paths[["dependency"]])
data.table::fwrite(joined, output_paths[["joined"]])

run_metadata <- data.frame(
  release = "unresolved_local_snapshot",
  model_metadata_file = basename(model_path),
  gene_effect_file = basename(effect_path),
  common_essential_file = basename(common_path),
  tnbc_models_measured = sum(model_manifest$model_group == "TNBC"),
  non_tnbc_breast_models_measured = sum(model_manifest$model_group == "Non_TNBC_breast"),
  dependency_threshold = config$depmap$dependency_threshold,
  model_probability_used = FALSE,
  stringsAsFactors = FALSE
)
data.table::fwrite(run_metadata, output_paths[["metadata"]])
invisible(write_checksum_manifest(
  unname(output_paths),
  file.path(analysis_root, "provenance", "depmap_outputs.sha256")
))

message("Measured TNBC models: ", run_metadata$tnbc_models_measured)
message("Measured non-TNBC breast models: ", run_metadata$non_tnbc_breast_models_measured)
message("DepMap genes summarized: ", nrow(dependency))
message("Primary DE rows with measured DepMap evidence: ", sum(!is.na(joined$tnbc_median_effect)))
