normalize_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\s*\\([0-9]+\\)$", "", x)
  x <- sub("_[0-9]+$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(x)
}

find_first_existing <- function(input_dir, candidates, required = FALSE) {
  paths <- file.path(input_dir, candidates)
  found <- paths[file.exists(paths)]
  if (length(found) > 0) return(found[[1]])
  if (required) stop("Missing required input; expected one of: ", paste(candidates, collapse = ", "))
  NA_character_
}

standardize_de_table <- function(de) {
  de <- suppressWarnings(janitor::clean_names(as.data.frame(de, check.names = FALSE)))
  gene_candidates <- c("gene_symbol", "symbol", "gene", "hugo_symbol", "gene_id")
  gene_col <- gene_candidates[gene_candidates %in% names(de)][1]
  lfc_candidates <- c("log2_fold_change", "log2fold_change", "log2foldchange", "log2_fc")
  lfc_col <- lfc_candidates[lfc_candidates %in% names(de)][1]
  if (is.na(gene_col)) stop("DE table has no recognized gene-symbol column")
  if (is.na(lfc_col) || !"padj" %in% names(de)) {
    stop("DE table must contain log2 fold-change and padj columns")
  }
  de$GENE <- normalize_gene_symbol(de[[gene_col]])
  de$LOG2FC <- as.numeric(de[[lfc_col]])
  de
}

filter_upregulated_de <- function(de, padj_threshold = 0.05, log2fc_threshold = 1) {
  de |>
    dplyr::filter(
      !is.na(.data$GENE), .data$GENE != "",
      !is.na(.data$padj), .data$padj < padj_threshold,
      !is.na(.data$LOG2FC), .data$LOG2FC > log2fc_threshold
    ) |>
    dplyr::arrange(.data$padj, dplyr::desc(.data$LOG2FC)) |>
    dplyr::distinct(.data$GENE, .keep_all = TRUE)
}

extract_model_vector <- function(matrix, model_id, value_name = "DEPENDENCY") {
  matrix <- as.data.frame(matrix, check.names = FALSE)
  raw_names <- names(matrix)
  cleaned_names <- names(suppressWarnings(janitor::clean_names(matrix)))
  id_candidates <- c("model_id", "modelid", "depmap_id", "v1", "x")
  id_clean <- id_candidates[id_candidates %in% cleaned_names][1]
  if (is.na(id_clean)) return(NULL)
  id_index <- match(id_clean, cleaned_names)
  selected <- matrix[tolower(matrix[[id_index]]) == tolower(model_id), , drop = FALSE]
  if (nrow(selected) == 0) return(NULL)
  gene_indices <- setdiff(seq_along(raw_names), id_index)
  gene_cols <- raw_names[gene_indices]
  values <- suppressWarnings(as.numeric(unlist(selected[1, gene_indices], use.names = FALSE)))
  result <- data.frame(
    GENE = normalize_gene_symbol(gene_cols),
    value = values,
    stringsAsFactors = FALSE
  )
  names(result)[2] <- value_name
  result <- result[!is.na(result[[value_name]]) & result$GENE != "", , drop = FALSE]
  result[!duplicated(result$GENE), , drop = FALSE]
}

join_de_with_model_effects <- function(de, effects, model_id) {
  de_up <- filter_upregulated_de(standardize_de_table(de))
  effect_vector <- extract_model_vector(effects, model_id, "DEPENDENCY")
  if (is.null(effect_vector)) stop("Model ID not found in gene-effect matrix: ", model_id)
  dplyr::inner_join(de_up, effect_vector, by = "GENE") |>
    dplyr::arrange(.data$DEPENDENCY)
}

prepare_plot_data <- function(joined, dependency_threshold = -0.5) {
  result <- joined |>
    dplyr::mutate(hit_base = .data$DEPENDENCY <= dependency_threshold)
  if ("P_DEP" %in% names(result)) {
    result <- result |>
      dplyr::mutate(
        dep50 = ifelse(
          !is.na(.data$P_DEP) & .data$P_DEP >= 0.5,
          "P(dep) >= 0.5", "P(dep) < 0.5 / NA"
        )
      )
  }
  result
}

find_model <- function(models, pattern = "mda.?mb.?231") {
  models <- suppressWarnings(janitor::clean_names(as.data.frame(models, check.names = FALSE)))
  name_candidates <- c("stripped_cell_line_name", "cell_line_name", "cell_line")
  id_candidates <- c("model_id", "modelid", "depmap_id")
  name_col <- name_candidates[name_candidates %in% names(models)][1]
  id_col <- id_candidates[id_candidates %in% names(models)][1]
  if (is.na(name_col) || is.na(id_col)) {
    stop("Model metadata lacks a recognized cell-line name or model ID column")
  }
  matched <- models[grepl(pattern, models[[name_col]], ignore.case = TRUE), , drop = FALSE]
  if (nrow(matched) == 0) stop("Could not find MDA-MB-231 in model metadata")
  list(id = as.character(matched[[id_col]][1]), name = as.character(matched[[name_col]][1]))
}

make_dependency_plot <- function(plot_data, output_file, dependency_threshold = -0.5) {
  color_column <- if ("dep50" %in% names(plot_data)) "dep50" else "hit_base"
  top_hits <- plot_data |>
    dplyr::filter(.data$DEPENDENCY <= dependency_threshold) |>
    dplyr::slice_min(.data$DEPENDENCY, n = 20)
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$LOG2FC, y = .data$DEPENDENCY, color = .data[[color_column]])
  ) +
    ggplot2::geom_hline(yintercept = dependency_threshold, linetype = 2, linewidth = 0.4) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "Overexpression vs Essentiality (MDA-MB-231)",
      subtitle = "DE: padj < 0.05 and log2FC > 1; candidate: Chronos <= -0.5",
      x = "log2 fold change (TNBC vs normal)",
      y = "CRISPR Chronos effect (lower = more essential)",
      color = NULL
    )
  if (nrow(top_hits) > 0) {
    plot <- plot + ggrepel::geom_text_repel(
      data = top_hits,
      ggplot2::aes(label = .data$GENE),
      size = 3, max.overlaps = Inf, seed = 42
    )
  }
  ggplot2::ggsave(output_file, plot, width = 8, height = 6, dpi = 300)
  invisible(plot)
}

run_depmap_integration <- function(
  input_dir = ".",
  output_dir = file.path(input_dir, "depmap_integration_outputs"),
  dependency_threshold = -0.5
) {
  de_file <- find_first_existing(input_dir, c("TNBC_vs_Normal_DEG.csv"), required = TRUE)
  effect_file <- find_first_existing(
    input_dir,
    c("CRISPRGeneEffect.csv", "CRISPR_gene_effect.csv", "Achilles_gene_effect.csv"),
    required = TRUE
  )
  model_file <- find_first_existing(
    input_dir,
    c("Model.csv", "Cell_lines.csv", "CellLines.csv", "sample_info.csv"),
    required = TRUE
  )
  probability_file <- find_first_existing(
    input_dir,
    c("CRISPRGeneDependency.csv", "CRISPR_gene_dependency.csv", "ScreenGeneDependency.csv")
  )
  common_essential_file <- find_first_existing(
    input_dir,
    c("CRISPRInferredCommonEssentials.csv", "AchillesCommonEssentialControls.csv")
  )
  de <- data.table::fread(de_file)
  effects <- data.table::fread(effect_file)
  models <- data.table::fread(model_file)
  model <- find_model(models)
  joined <- join_de_with_model_effects(de, effects, model$id)
  if (!is.na(probability_file)) {
    probabilities <- data.table::fread(probability_file)
    probability_vector <- extract_model_vector(probabilities, model$id, "P_DEP")
    if (is.null(probability_vector)) {
      message("Probability file is screen-level or lacks model ", model$id, "; continuing without P_DEP.")
    } else {
      joined <- dplyr::left_join(joined, probability_vector, by = "GENE")
    }
  }
  if (!is.na(common_essential_file)) {
    common_essential_table <- data.table::fread(common_essential_file, sep = ",")
    common_essential_genes <- normalize_gene_symbol(common_essential_table[[1]])
    joined$COMMON_ESSENTIAL <- joined$GENE %in% common_essential_genes
  } else {
    joined$COMMON_ESSENTIAL <- NA
  }
  candidates <- joined |>
    dplyr::filter(.data$DEPENDENCY <= dependency_threshold) |>
    dplyr::arrange(.data$DEPENDENCY, dplyr::desc(.data$LOG2FC))
  noncommon_candidates <- candidates |>
    dplyr::filter(!is.na(.data$COMMON_ESSENTIAL), !.data$COMMON_ESSENTIAL)
  plot_data <- prepare_plot_data(joined, dependency_threshold)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  joined_file <- file.path(output_dir, "DE_up_vs_DepMap_MDA_MB_231_joined.csv")
  candidate_file <- file.path(output_dir, "Candidates_overexpressed_and_essential_MDA_MB_231.csv")
  noncommon_candidate_file <- file.path(
    output_dir,
    "Candidates_noncommon_overexpressed_and_essential_MDA_MB_231.csv"
  )
  plot_file <- file.path(output_dir, "Scatter_log2FC_vs_Dependency_MDA_MB_231.png")
  data.table::fwrite(joined, joined_file)
  data.table::fwrite(candidates, candidate_file)
  data.table::fwrite(noncommon_candidates, noncommon_candidate_file)
  make_dependency_plot(plot_data, plot_file, dependency_threshold)
  message("Found MDA-MB-231: ", model$name, " | Model ID: ", model$id)
  message("Matched DE/DepMap genes: ", nrow(joined))
  message("Candidates at Chronos <= ", dependency_threshold, ": ", nrow(candidates))
  message("Candidates not marked common essential: ", nrow(noncommon_candidates))
  list(
    model = model,
    joined = joined,
    candidates = candidates,
    noncommon_candidates = noncommon_candidates,
    files = c(joined_file, candidate_file, noncommon_candidate_file, plot_file)
  )
}
