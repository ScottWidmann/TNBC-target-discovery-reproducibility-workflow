#!/usr/bin/env Rscript

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_argument) != 1) {
  stop("Run this stage with Rscript", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
analysis_root <- file.path(project_root, "analysis_v2")
results_dir <- file.path(analysis_root, "results")
plots_dir <- file.path(analysis_root, "plots")
provenance_dir <- file.path(analysis_root, "provenance")

source(file.path(analysis_root, "R", "io_contracts.R"))
source(file.path(analysis_root, "R", "evidence.R"))
source(file.path(analysis_root, "R", "reporting.R"))

config <- read_analysis_config(file.path(analysis_root, "config", "analysis.yml"))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(name) {
  path <- file.path(results_dir, name)
  if (!file.exists(path)) {
    stop("Required Stage 05 input is missing: ", path, call. = FALSE)
  }
  data.table::fread(path, data.table = FALSE)
}

primary <- read_result("de_primary_tnbc_vs_normal.csv")
paired <- read_result("de_paired_tnbc_vs_normal.csv")
specificity <- read_result("de_specificity_tnbc_vs_other.csv")
depmap <- read_result("depmap_gene_dependency_summary.csv")
tractability <- data.table::fread(
  file.path(analysis_root, "config", "tractability_evidence.csv"),
  data.table = FALSE
)

candidate_evidence <- assemble_candidate_evidence(
  primary,
  paired,
  specificity,
  depmap,
  tractability,
  config
)
candidate_path <- file.path(results_dir, "candidate_evidence.csv")
data.table::fwrite(candidate_evidence, candidate_path)

manifest <- read_result("sample_manifest_included.csv")
counts_path <- file.path(results_dir, "intermediate", "counts.rds")
if (!file.exists(counts_path)) {
  stop("Required count matrix is missing: ", counts_path, call. = FALSE)
}
counts <- readRDS(counts_path)
primary_manifest <- manifest[
  manifest$analysis_group %in% c("TNBC", "Normal"),
  ,
  drop = FALSE
]
primary_counts <- counts[, primary_manifest$sample_id, drop = FALSE]
set.seed(config$random_seed)
pca <- compute_pca_scores(primary_counts, primary_manifest, n_top = 500)
group_counts <- table(primary_manifest$analysis_group)

plot_paths <- c(
  pca = file.path(plots_dir, "pca_primary.png"),
  volcano = file.path(plots_dir, "volcano_primary.png"),
  paired = file.path(plots_dir, "concordance_paired_primary.png"),
  specificity = file.path(plots_dir, "concordance_specificity_primary.png"),
  depmap = file.path(plots_dir, "depmap_tnbc_vs_non_tnbc.png")
)
invisible(save_report_plot(make_pca_plot(pca, group_counts), plot_paths[["pca"]], 8, 6))
invisible(save_report_plot(
  make_volcano_plot(
    primary,
    config$differential_expression$alpha,
    config$differential_expression$lfc_threshold
  ),
  plot_paths[["volcano"]],
  8,
  6
))
invisible(save_report_plot(
  make_concordance_plot(primary, paired, "Paired"),
  plot_paths[["paired"]],
  8,
  6
))
invisible(save_report_plot(
  make_concordance_plot(primary, specificity, "Specificity"),
  plot_paths[["specificity"]],
  8,
  6
))

depmap_metadata <- read_result("depmap_run_metadata.csv")
invisible(save_report_plot(
  make_depmap_plot(
    depmap,
    candidate_evidence,
    config,
    depmap_metadata$tnbc_models_measured[[1]],
    depmap_metadata$non_tnbc_breast_models_measured[[1]]
  ),
  plot_paths[["depmap"]],
  8,
  6
))

enrichment_inputs <- c(
  enrichment_go_up = "ora_go_bp_up.csv",
  enrichment_go_down = "ora_go_bp_down.csv",
  enrichment_kegg_up = "ora_kegg_legacy_up.csv",
  enrichment_kegg_down = "ora_kegg_legacy_down.csv"
)
enrichment_titles <- c(
  enrichment_go_up = "GO biological process: upregulated genes",
  enrichment_go_down = "GO biological process: downregulated genes",
  enrichment_kegg_up = "KEGG legacy: upregulated genes",
  enrichment_kegg_down = "KEGG legacy: downregulated genes"
)
for (plot_name in names(enrichment_inputs)) {
  enrichment <- read_result(enrichment_inputs[[plot_name]])
  enrichment_plot <- make_enrichment_dotplot(
    enrichment,
    enrichment_titles[[plot_name]],
    n_tnbc = unique(primary$n_numerator)[[1]],
    n_normal = unique(primary$n_denominator)[[1]],
    alpha = config$differential_expression$alpha,
    lfc_threshold = config$differential_expression$lfc_threshold
  )
  if (!is.null(enrichment_plot)) {
    path <- file.path(plots_dir, paste0(plot_name, ".png"))
    invisible(save_report_plot(enrichment_plot, path, 10, 7.5))
    plot_paths[[plot_name]] <- path
  }
}

warnings <- data.frame(
  warning_id = sprintf("W%03d", 1:7),
  scope = c(
    "DepMap provenance", "Tractability", "Scientific claim", "Cohort",
    "Tumor subtype", "Model adjustment", "Identifier mapping"
  ),
  warning = c(
    "The exact DepMap release is unresolved for the frozen local snapshot.",
    "Tractability is not assessed because no versioned evidence source was supplied.",
    "This observational workflow nominates targets; it does not establish causal dependency or therapeutic efficacy.",
    "The primary analysis uses all eligible normals; only 11 TNBC-normal patient pairs support the paired sensitivity analysis.",
    "TNBC is defined from clinical ER, PR, and HER2 fields, not a molecular-subtype classifier.",
    "Age was eligible for primary and specificity adjustment; race and tissue source site were rejected as group-confounded.",
    "Enrichment is limited to tested genes that map to the frozen local GO and KEGG legacy collections."
  ),
  stringsAsFactors = FALSE
)
warnings_path <- file.path(provenance_dir, "warnings.csv")
data.table::fwrite(warnings, warnings_path)

session_path <- file.path(provenance_dir, "session_info.txt")
writeLines(capture.output(utils::sessionInfo()), session_path, useBytes = TRUE)

source_paths <- c(
  list.files(file.path(analysis_root, "R"), recursive = TRUE, full.names = TRUE),
  list.files(file.path(analysis_root, "scripts"), recursive = TRUE, full.names = TRUE),
  list.files(file.path(analysis_root, "config"), recursive = TRUE, full.names = TRUE),
  list.files(file.path(analysis_root, "tests"), recursive = TRUE, full.names = TRUE),
  file.path(analysis_root, "run_all.R"),
  file.path(analysis_root, "corrected_analysis_report.Rmd"),
  file.path(project_root, "docs", "superpowers", "specs", "2026-08-19-corrected-tnbc-analysis-design.md"),
  file.path(project_root, "docs", "superpowers", "plans", "2026-08-19-corrected-tnbc-analysis.md")
)
source_paths <- sort(unique(source_paths[file.exists(source_paths)]))
invisible(write_checksum_manifest(
  source_paths,
  file.path(provenance_dir, "source_checksums.sha256")
))

report_path <- file.path(results_dir, "corrected_analysis_report.html")
invisible(render_corrected_report(report_path))
invisible(validate_rendered_report(report_path, minimum_images = 9L))

claimed_outputs <- c(
  list.files(results_dir, recursive = TRUE, full.names = TRUE),
  list.files(plots_dir, recursive = TRUE, full.names = TRUE),
  warnings_path,
  session_path,
  file.path(provenance_dir, "source_checksums.sha256")
)
claimed_outputs <- setdiff(
  claimed_outputs,
  file.path(provenance_dir, "output_manifest.sha256")
)
invisible(write_output_manifest(
  claimed_outputs,
  file.path(provenance_dir, "output_manifest.sha256")
))
invisible(validate_final_outputs(analysis_root))

class_counts <- table(candidate_evidence$evidence_class)
message("Candidate evidence rows: ", nrow(candidate_evidence))
message(
  "Evidence classes: ",
  paste(names(class_counts), as.integer(class_counts), sep = "=", collapse = ", ")
)
message("Rendered report: ", report_path)
