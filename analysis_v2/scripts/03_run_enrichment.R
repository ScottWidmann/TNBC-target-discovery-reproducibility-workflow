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
source(file.path(analysis_root, "R", "enrichment.R"))

config <- read_analysis_config(file.path(analysis_root, "config", "analysis.yml"))
results_dir <- file.path(analysis_root, "results")
de <- data.table::fread(
  file.path(results_dir, "de_primary_tnbc_vs_normal.csv"),
  data.table = FALSE
)
de_config <- config$differential_expression

message("Building local GO Biological Process mapping")
go_mapping <- build_go_bp_mapping()
message("Building local MSigDB KEGG-legacy mapping")
kegg_mapping <- build_kegg_legacy_mapping()

databases <- list(GO_BP = go_mapping, KEGG_LEGACY = kegg_mapping)
ora_results <- list()
ranked_results <- list()
audit_rows <- list()

for (database in names(databases)) {
  mapping <- databases[[database]]
  mapping_audit <- map_tested_genes(de, mapping, database)$audit
  mapping_audit$database_version <- paste(unique(mapping$database_version), collapse = ";")
  audit_rows[[paste0(database, "_mapping")]] <- mapping_audit

  for (direction in c("up", "down")) {
    result <- run_ora(
      de,
      mapping,
      direction = direction,
      alpha = de_config$alpha,
      lfc_threshold = de_config$lfc_threshold
    )
    result$database <- database
    ora_results[[paste(database, direction, sep = "_")]] <- result
    audit_rows[[paste(database, direction, sep = "_")]] <- data.frame(
      database = database,
      direction = direction,
      status = attr(result, "status"),
      foreground_size = attr(result, "foreground_size"),
      universe_size = attr(result, "universe_size"),
      enriched_terms_fdr_0_05 = sum(result$padj < 0.05, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  ranked <- run_ranked_enrichment(
    de,
    mapping,
    seed = config$random_seed
  )
  ranked$database <- database
  ranked_results[[database]] <- ranked
}

output_paths <- character()
for (name in names(ora_results)) {
  path <- file.path(results_dir, paste0("ora_", tolower(name), ".csv"))
  data.table::fwrite(ora_results[[name]], path)
  output_paths <- c(output_paths, path)
}
for (name in names(ranked_results)) {
  path <- file.path(results_dir, paste0("gsea_", tolower(name), ".csv"))
  data.table::fwrite(ranked_results[[name]], path)
  output_paths <- c(output_paths, path)
}

mapping_audits <- data.table::rbindlist(audit_rows, fill = TRUE)
mapping_audit_path <- file.path(results_dir, "mapping_audit.csv")
data.table::fwrite(mapping_audits, mapping_audit_path)
output_paths <- c(output_paths, mapping_audit_path)

metadata <- data.frame(
  analysis = c("GO_BP", "KEGG_LEGACY"),
  source = c("org.Hs.eg.db plus GO.db", "msigdbr C2:CP:KEGG_LEGACY"),
  source_version = c(
    paste0(
      "org.Hs.eg.db ", utils::packageVersion("org.Hs.eg.db"),
      "; GO.db ", utils::packageVersion("GO.db")
    ),
    paste(unique(kegg_mapping$database_version), collapse = ";")
  ),
  online_lookup = FALSE,
  random_seed = config$random_seed,
  stringsAsFactors = FALSE
)
metadata_path <- file.path(results_dir, "enrichment_run_metadata.csv")
data.table::fwrite(metadata, metadata_path)
output_paths <- c(output_paths, metadata_path)

invisible(write_checksum_manifest(
  output_paths,
  file.path(analysis_root, "provenance", "enrichment_outputs.sha256")
))
print(mapping_audits)
