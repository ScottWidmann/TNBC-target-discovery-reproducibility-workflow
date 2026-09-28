normalize_depmap_gene_symbol <- function(x) {
  value <- trimws(as.character(x))
  value <- sub("\\s*\\([0-9]+\\)$", "", value)
  toupper(value)
}

classify_depmap_models <- function(models) {
  require_columns(
    models,
    c("ModelID", "CellLineName", "OncotreeLineage", "ModelSubtypeFeatures"),
    "DepMap model metadata"
  )
  lineage <- toupper(trimws(as.character(models$OncotreeLineage)))
  features <- toupper(trimws(as.character(models$ModelSubtypeFeatures)))
  features[is.na(features)] <- ""
  is_breast <- lineage == "BREAST"
  is_tnbc <- is_breast & grepl("\\bTNBC\\b", features)
  is_non_tnbc <- is_breast & !is_tnbc & grepl("ER\\+|HER2\\+", features)
  keep <- is_tnbc | is_non_tnbc
  data.frame(
    model_id = as.character(models$ModelID[keep]),
    model_name = as.character(models$CellLineName[keep]),
    model_group = ifelse(is_tnbc[keep], "TNBC", "Non_TNBC_breast"),
    source_label = as.character(models$ModelSubtypeFeatures[keep]),
    stringsAsFactors = FALSE
  )
}

depmap_model_id_column <- function(effects) {
  candidates <- c("ModelID", "model_id", "DepMap_ID", "depmap_id", "V1", "X")
  selected <- candidates[candidates %in% names(effects)]
  if (length(selected) == 0) {
    stop("DepMap gene-effect matrix has no recognized model ID column", call. = FALSE)
  }
  selected[[1]]
}

extract_depmap_effects <- function(effects, model_manifest) {
  require_columns(
    model_manifest,
    c("model_id", "model_name", "model_group"),
    "DepMap model manifest"
  )
  id_column <- depmap_model_id_column(effects)
  effect_model_ids <- as.character(effects[[id_column]])
  missing_models <- setdiff(model_manifest$model_id, effect_model_ids)
  if (length(missing_models) > 0) {
    stop(
      "DepMap gene-effect matrix is missing panel models: ",
      paste(missing_models, collapse = ", "),
      call. = FALSE
    )
  }
  gene_columns <- setdiff(names(effects), id_column)
  gene_symbols <- normalize_depmap_gene_symbol(gene_columns)
  duplicate_symbols <- unique(gene_symbols[duplicated(gene_symbols)])
  if (length(duplicate_symbols) > 0) {
    stop(
      "DepMap gene-effect columns collapse to duplicate symbols: ",
      paste(duplicate_symbols, collapse = ", "),
      call. = FALSE
    )
  }
  selected <- effects[
    match(model_manifest$model_id, effect_model_ids),
    c(id_column, gene_columns),
    drop = FALSE
  ]
  names(selected)[names(selected) == id_column] <- "model_id"
  selected <- data.table::as.data.table(selected)
  long <- data.table::melt(
    selected,
    id.vars = "model_id",
    measure.vars = gene_columns,
    variable.name = "gene_column",
    value.name = "effect",
    variable.factor = FALSE
  )
  gene_lookup <- stats::setNames(gene_symbols, gene_columns)
  long[, gene_symbol := unname(gene_lookup[gene_column])]
  long[, model_name := model_manifest$model_name[match(model_id, model_manifest$model_id)]]
  long[, model_group := model_manifest$model_group[match(model_id, model_manifest$model_id)]]
  long[, effect := suppressWarnings(as.numeric(effect))]
  long <- long[is.finite(effect), .(
    model_id = as.character(model_id),
    model_name = as.character(model_name),
    model_group = as.character(model_group),
    gene_symbol = as.character(gene_symbol),
    effect
  )]
  data.table::setorder(long, model_id, gene_symbol)
  as.data.frame(long)
}

summarize_depmap_by_gene <- function(effects, dependency_threshold = -0.5) {
  require_columns(
    effects,
    c("gene_symbol", "model_name", "model_group", "effect"),
    "long DepMap effects"
  )
  summarize_group <- function(group_name, prefix) {
    selected <- effects[
      effects$model_group == group_name & is.finite(effects$effect),
      ,
      drop = FALSE
    ]
    summary <- data.table::as.data.table(selected)[, .(
      n = as.integer(.N),
      median_effect = stats::median(effect),
      q1_effect = as.numeric(stats::quantile(effect, 0.25, names = FALSE)),
      q3_effect = as.numeric(stats::quantile(effect, 0.75, names = FALSE)),
      fraction_below_threshold = mean(effect <= dependency_threshold)
    ), by = gene_symbol]
    data.table::setnames(
      summary,
      c("n", "median_effect", "q1_effect", "q3_effect", "fraction_below_threshold"),
      paste0(prefix, c(
        "_n", "_median_effect", "_q1_effect", "_q3_effect",
        "_fraction_below_threshold"
      ))
    )
    summary
  }

  tnbc <- summarize_group("TNBC", "tnbc")
  non_tnbc <- summarize_group("Non_TNBC_breast", "non_tnbc")
  summary <- merge(tnbc, non_tnbc, by = "gene_symbol", all = TRUE, sort = TRUE)
  summary$tnbc_minus_non_tnbc <-
    summary$tnbc_median_effect - summary$non_tnbc_median_effect

  normalized_names <- toupper(gsub("[^A-Z0-9]", "", effects$model_name))
  mda <- effects[normalized_names == "MDAMB231", c("gene_symbol", "effect"), drop = FALSE]
  if (nrow(mda) > 0) {
    if (anyDuplicated(mda$gene_symbol)) {
      stop("MDA-MB-231 has duplicated gene-effect rows", call. = FALSE)
    }
    summary$mda_mb_231_effect <- mda$effect[match(summary$gene_symbol, mda$gene_symbol)]
  } else {
    summary$mda_mb_231_effect <- NA_real_
  }
  as.data.frame(summary)
}

annotate_common_essentials <- function(summary, common_table) {
  require_columns(summary, "gene_symbol", "DepMap dependency summary")
  result <- summary
  if (is.null(common_table)) {
    result$common_essential <- NA
    return(result)
  }
  if (ncol(common_table) == 0) {
    stop("Common-essential table has no columns", call. = FALSE)
  }
  common_genes <- normalize_depmap_gene_symbol(common_table[[1]])
  common_genes <- unique(common_genes[!is.na(common_genes) & common_genes != ""])
  result$common_essential <- normalize_depmap_gene_symbol(result$gene_symbol) %in% common_genes
  result
}

read_common_essential_table <- function(path) {
  if (!file.exists(path)) {
    stop("Common-essential file does not exist: ", path, call. = FALSE)
  }
  table <- data.table::fread(
    path,
    sep = ",",
    header = TRUE,
    data.table = FALSE
  )
  if (ncol(table) != 1) {
    stop("Common-essential file must contain exactly one column", call. = FALSE)
  }
  table
}

join_depmap_to_de <- function(de, dependency) {
  require_columns(de, c("ensembl_id", "gene_symbol"), "differential-expression results")
  require_columns(dependency, "gene_symbol", "DepMap dependency summary")
  if (anyDuplicated(de$ensembl_id)) {
    stop("Differential-expression results contain duplicate Ensembl IDs", call. = FALSE)
  }
  normalized_dependency <- normalize_depmap_gene_symbol(dependency$gene_symbol)
  if (anyDuplicated(normalized_dependency)) {
    stop("DepMap dependency summary contains duplicate gene symbols", call. = FALSE)
  }
  dependency_index <- match(
    normalize_depmap_gene_symbol(de$gene_symbol),
    normalized_dependency
  )
  evidence_columns <- setdiff(names(dependency), "gene_symbol")
  evidence <- dependency[dependency_index, evidence_columns, drop = FALSE]
  rownames(evidence) <- NULL
  data.frame(de, evidence, stringsAsFactors = FALSE, check.names = FALSE)
}
