assert_unique_key <- function(x, key, object_name) {
  require_columns(x, key, object_name)
  values <- as.character(x[[key]])
  if (anyNA(values) || any(!nzchar(values))) {
    stop(object_name, " contains a missing ", key, call. = FALSE)
  }
  if (anyDuplicated(values)) {
    stop(object_name, " contains duplicate ", key, " values", call. = FALSE)
  }
  invisible(TRUE)
}

logical_true <- function(x) {
  !is.na(x) & x
}

assign_evidence_class <- function(x, config) {
  required <- c(
    "primary_hit", "paired_concordant", "specificity_concordant",
    "tnbc_dependency", "depmap_selective", "common_essential"
  )
  require_columns(x, required, "candidate evidence")

  primary <- logical_true(x$primary_hit)
  paired <- logical_true(x$paired_concordant)
  specificity <- logical_true(x$specificity_concordant)
  dependency <- logical_true(x$tnbc_dependency)
  selective <- logical_true(x$depmap_selective)
  not_common <- !is.na(x$common_essential) & !x$common_essential

  class <- rep("not_candidate", nrow(x))
  class[primary] <- "C_expression_only"
  class[primary & dependency & (paired | specificity)] <-
    "B_expression_dependency"
  class[
    primary & paired & specificity & dependency & selective & not_common
  ] <- "A_multistream"
  class
}

join_sensitivity_fields <- function(out, sensitivity, prefix) {
  assert_unique_key(sensitivity, "ensembl_id", paste(prefix, "DE results"))
  require_columns(
    sensitivity,
    c("padj", "log2fc_shrunk"),
    paste(prefix, "DE results")
  )
  index <- match(out$ensembl_id, sensitivity$ensembl_id)
  out[[paste0(prefix, "_padj")]] <- sensitivity$padj[index]
  out[[paste0(prefix, "_log2fc_shrunk")]] <- sensitivity$log2fc_shrunk[index]
  out
}

join_depmap_fields <- function(out, depmap) {
  required <- c(
    "gene_symbol", "tnbc_median_effect",
    "tnbc_fraction_below_threshold", "tnbc_minus_non_tnbc",
    "common_essential"
  )
  require_columns(depmap, required, "DepMap summary")
  if (anyDuplicated(depmap$gene_symbol)) {
    stop("DepMap summary contains duplicate gene_symbol values", call. = FALSE)
  }

  index <- match(toupper(out$gene_symbol), toupper(depmap$gene_symbol))
  fields <- setdiff(names(depmap), "gene_symbol")
  for (field in fields) {
    out[[field]] <- depmap[[field]][index]
  }
  out
}

join_tractability_fields <- function(out, tractability) {
  required <- c(
    "ensembl_id", "gene_symbol", "tractability_status", "source_name",
    "source_version", "source_record_id", "evidence_note"
  )
  require_columns(tractability, required, "tractability evidence")

  fields <- setdiff(required, c("ensembl_id", "gene_symbol"))
  if (!nrow(tractability)) {
    for (field in fields) {
      out[[paste0("tractability_", field)]] <- NA_character_
    }
    out$tractability_assessment <- "not_assessed"
    return(out)
  }

  if (anyDuplicated(tractability$ensembl_id)) {
    stop("tractability evidence contains duplicate ensembl_id values", call. = FALSE)
  }
  index <- match(out$ensembl_id, tractability$ensembl_id)
  for (field in fields) {
    out[[paste0("tractability_", field)]] <- tractability[[field]][index]
  }
  out$tractability_assessment <- ifelse(
    is.na(out$tractability_tractability_status) |
      !nzchar(out$tractability_tractability_status),
    "not_assessed",
    "assessed"
  )
  out
}

assemble_candidate_evidence <- function(
    primary,
    paired,
    specificity,
    depmap,
    tractability,
    config) {
  assert_unique_key(primary, "ensembl_id", "primary DE results")
  require_columns(
    primary,
    c("gene_symbol", "gene_type", "padj", "log2fc_shrunk"),
    "primary DE results"
  )

  identity_fields <- intersect(
    c("ensembl_id", "gene_symbol", "gene_type"),
    names(primary)
  )
  out <- primary[, identity_fields, drop = FALSE]
  primary_fields <- setdiff(names(primary), identity_fields)
  for (field in primary_fields) {
    out[[paste0("primary_", field)]] <- primary[[field]]
  }

  out <- join_sensitivity_fields(out, paired, "paired")
  out <- join_sensitivity_fields(out, specificity, "specificity")
  out <- join_depmap_fields(out, depmap)
  out <- join_tractability_fields(out, tractability)

  alpha <- config$differential_expression$alpha
  lfc_threshold <- config$differential_expression$lfc_threshold
  dependency_threshold <- config$depmap$dependency_threshold
  dependency_fraction <- config$depmap$dependency_fraction
  selectivity_delta <- config$depmap$selectivity_delta

  out$primary_hit <- !is.na(out$primary_padj) &
    out$primary_padj < alpha &
    !is.na(out$primary_log2fc_shrunk) &
    out$primary_log2fc_shrunk > lfc_threshold
  out$paired_available <- !is.na(out$paired_log2fc_shrunk)
  out$paired_concordant <- out$paired_available & out$paired_log2fc_shrunk > 0
  out$specificity_available <- !is.na(out$specificity_log2fc_shrunk)
  out$specificity_concordant <- out$specificity_available &
    out$specificity_log2fc_shrunk > 0
  out$depmap_available <- !is.na(out$tnbc_median_effect) &
    !is.na(out$tnbc_fraction_below_threshold)
  out$tnbc_dependency <- out$depmap_available &
    out$tnbc_median_effect <= dependency_threshold &
    out$tnbc_fraction_below_threshold >= dependency_fraction
  out$depmap_selectivity_available <- !is.na(out$tnbc_minus_non_tnbc)
  out$depmap_selective <- out$depmap_selectivity_available &
    out$tnbc_minus_non_tnbc <= selectivity_delta
  out$evidence_class <- assign_evidence_class(out, config)

  class_order <- c(
    "A_multistream", "B_expression_dependency",
    "C_expression_only", "not_candidate"
  )
  out$evidence_class <- factor(out$evidence_class, levels = class_order)
  order_index <- order(
    out$evidence_class,
    out$primary_padj,
    -out$primary_log2fc_shrunk,
    out$ensembl_id,
    na.last = TRUE
  )
  out <- out[order_index, , drop = FALSE]
  out$evidence_class <- as.character(out$evidence_class)
  rownames(out) <- NULL
  out
}
