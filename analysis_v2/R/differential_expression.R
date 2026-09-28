filter_counts_for_contrast <- function(
  counts,
  sample_ids,
  min_count,
  min_fraction
) {
  missing_samples <- setdiff(sample_ids, colnames(counts))
  if (length(missing_samples) > 0) {
    stop(
      "Count matrix is missing requested samples: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(min_count) != 1 || !is.finite(min_count) || min_count < 0) {
    stop("min_count must be one nonnegative number", call. = FALSE)
  }
  if (length(min_fraction) != 1 || !is.finite(min_fraction) ||
      min_fraction <= 0 || min_fraction > 1) {
    stop("min_fraction must be greater than 0 and at most 1", call. = FALSE)
  }
  selected <- counts[, sample_ids, drop = FALSE]
  minimum_samples <- ceiling(min_fraction * length(sample_ids))
  keep <- rowSums(selected >= min_count) >= minimum_samples
  selected[keep, , drop = FALSE]
}

join_gene_annotation <- function(de, genes) {
  require_columns(de, "ensembl_id", "differential-expression results")
  require_columns(
    genes,
    c("ensembl_id", "gene_symbol", "gene_type"),
    "gene annotation"
  )
  if (anyDuplicated(de$ensembl_id)) {
    stop("Differential-expression results contain duplicate Ensembl IDs", call. = FALSE)
  }
  if (anyDuplicated(genes$ensembl_id)) {
    stop("Gene annotation contains duplicate Ensembl IDs", call. = FALSE)
  }
  gene_index <- match(de$ensembl_id, genes$ensembl_id)
  if (anyNA(gene_index)) {
    stop("Gene annotation is missing tested Ensembl IDs", call. = FALSE)
  }
  annotation <- genes[gene_index, c("gene_symbol", "gene_type"), drop = FALSE]
  data.frame(
    ensembl_id = de$ensembl_id,
    annotation,
    de[setdiff(names(de), "ensembl_id")],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

nonmissing_covariate <- function(x) {
  present <- !is.na(x)
  if (is.character(x) || is.factor(x)) {
    value <- trimws(as.character(x))
    present <- present & value != "" &
      !toupper(value) %in% c("[NOT AVAILABLE]", "[UNKNOWN]", "UNKNOWN")
  }
  present
}

assess_covariates <- function(col_data, candidates, min_complete = 0.9) {
  require_columns(col_data, "group", "DESeq2 column data")
  groups <- unique(as.character(col_data$group))
  if (length(groups) != 2) {
    stop("Covariate assessment requires exactly two groups", call. = FALSE)
  }
  audits <- lapply(candidates, function(candidate) {
    if (!candidate %in% names(col_data)) {
      return(data.frame(
        covariate = candidate,
        completeness = 0,
        eligible = FALSE,
        reason = "absent_from_column_data"
      ))
    }
    value <- col_data[[candidate]]
    present <- nonmissing_covariate(value)
    completeness <- mean(present)
    if (completeness < min_complete) {
      return(data.frame(
        covariate = candidate,
        completeness = completeness,
        eligible = FALSE,
        reason = "below_completeness_threshold"
      ))
    }
    complete_data <- col_data[present, , drop = FALSE]
    complete_value <- value[present]
    if (is.numeric(complete_value) || is.integer(complete_value)) {
      variation <- vapply(
        split(complete_value, as.character(complete_data$group)),
        function(x) length(unique(x)) >= 2,
        FUN.VALUE = logical(1)
      )
      if (!all(variation)) {
        return(data.frame(
          covariate = candidate,
          completeness = completeness,
          eligible = FALSE,
          reason = "no_within_group_variation"
        ))
      }
    } else {
      complete_value <- droplevels(factor(complete_value))
      if (nlevels(complete_value) < 2) {
        return(data.frame(
          covariate = candidate,
          completeness = completeness,
          eligible = FALSE,
          reason = "fewer_than_two_levels"
        ))
      }
      cross_tabulation <- table(complete_data$group, complete_value)
      if (any(colSums(cross_tabulation > 0) < 2)) {
        return(data.frame(
          covariate = candidate,
          completeness = completeness,
          eligible = FALSE,
          reason = "confounded_with_group"
        ))
      }
      if (any(colSums(cross_tabulation) < 5)) {
        return(data.frame(
          covariate = candidate,
          completeness = completeness,
          eligible = FALSE,
          reason = "sparse_covariate_level"
        ))
      }
    }
    design_data <- data.frame(
      covariate = complete_value,
      group = factor(complete_data$group)
    )
    design <- stats::model.matrix(~ covariate + group, design_data)
    if (qr(design)$rank != ncol(design)) {
      return(data.frame(
        covariate = candidate,
        completeness = completeness,
        eligible = FALSE,
        reason = "rank_deficient_with_group"
      ))
    }
    data.frame(
      covariate = candidate,
      completeness = completeness,
      eligible = TRUE,
      reason = "eligible"
    )
  })
  do.call(rbind, audits)
}

build_deseq_design <- function(covariates = character(), paired = FALSE) {
  if (paired && length(covariates) > 0) {
    stop("Paired models do not accept additional covariates", call. = FALSE)
  }
  if (paired) {
    return(stats::as.formula("~ patient_id + group"))
  }
  terms <- c(covariates, "group")
  stats::as.formula(paste("~", paste(terms, collapse = " + ")))
}

prepare_model_covariates <- function(col_data, covariates) {
  result <- col_data
  for (covariate in covariates) {
    value <- result[[covariate]]
    if (is.numeric(value) || is.integer(value)) {
      standard_deviation <- stats::sd(value)
      if (!is.finite(standard_deviation) || standard_deviation == 0) {
        stop(
          "Numeric covariate has no usable variation: ", covariate,
          call. = FALSE
        )
      }
      result[[covariate]] <- as.numeric(scale(value))
    } else {
      result[[covariate]] <- droplevels(factor(value))
    }
  }
  result
}

encode_patient_block <- function(patient_id) {
  patient_id <- as.character(patient_id)
  codes <- match(patient_id, unique(patient_id))
  labels <- sprintf("patient_%04d", codes)
  factor(labels, levels = unique(labels))
}

run_deseq_contrast <- function(
  counts,
  col_data,
  numerator,
  denominator,
  covariates = character(),
  paired = FALSE,
  shrink_type = "apeglm"
) {
  required_columns <- c("sample_id", "group")
  if (paired) required_columns <- c(required_columns, "patient_id")
  require_columns(col_data, required_columns, "DESeq2 column data")
  requested <- as.character(col_data$group) %in% c(numerator, denominator)
  model_data <- col_data[requested, , drop = FALSE]
  if (nrow(model_data) == 0) {
    stop("No samples found for requested contrast", call. = FALSE)
  }
  missing_count_samples <- setdiff(model_data$sample_id, colnames(counts))
  if (length(missing_count_samples) > 0) {
    stop(
      "Count matrix is missing DESeq2 samples: ",
      paste(missing_count_samples, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(covariates) > 0) {
    require_columns(model_data, covariates, "DESeq2 column data")
    complete <- stats::complete.cases(model_data[, covariates, drop = FALSE])
    model_data <- model_data[complete, , drop = FALSE]
    model_data <- prepare_model_covariates(model_data, covariates)
  }
  group_counts <- table(model_data$group)
  if (!all(c(numerator, denominator) %in% names(group_counts)) ||
      any(group_counts[c(numerator, denominator)] < 2)) {
    stop("Each contrast group must contain at least two samples", call. = FALSE)
  }
  if (paired) {
    pair_table <- table(model_data$patient_id, model_data$group)
    required_columns_in_pair <- match(c(denominator, numerator), colnames(pair_table))
    complete_pairs <- rowSums(pair_table[, required_columns_in_pair, drop = FALSE] == 1) == 2
    model_data <- model_data[model_data$patient_id %in% rownames(pair_table)[complete_pairs], , drop = FALSE]
    if (length(unique(model_data$patient_id)) < 2) {
      stop("Paired contrast requires at least two complete patient pairs", call. = FALSE)
    }
  }
  model_data$group <- factor(model_data$group, levels = c(denominator, numerator))
  if (paired) model_data$patient_id <- encode_patient_block(model_data$patient_id)
  rownames(model_data) <- model_data$sample_id
  model_counts <- counts[, model_data$sample_id, drop = FALSE]
  design <- build_deseq_design(covariates, paired)
  design_matrix <- stats::model.matrix(design, model_data)
  if (qr(design_matrix)$rank != ncol(design_matrix)) {
    stop("DESeq2 design matrix is not full rank", call. = FALSE)
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = model_counts,
    colData = model_data,
    design = design
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  raw <- DESeq2::results(dds, contrast = c("group", numerator, denominator))
  expected_coefficient <- paste0("group_", numerator, "_vs_", denominator)
  coefficients <- DESeq2::resultsNames(dds)
  if (!expected_coefficient %in% coefficients) {
    stop(
      "Expected DESeq2 coefficient not found: ", expected_coefficient,
      "; available: ", paste(coefficients, collapse = ", "),
      call. = FALSE
    )
  }
  shrunken <- DESeq2::lfcShrink(
    dds,
    coef = expected_coefficient,
    type = shrink_type,
    quiet = TRUE
  )
  data.frame(
    ensembl_id = rownames(raw),
    base_mean = raw$baseMean,
    log2fc_raw = raw$log2FoldChange,
    lfc_se_raw = raw$lfcSE,
    stat = raw$stat,
    pvalue = raw$pvalue,
    padj = raw$padj,
    log2fc_shrunk = shrunken$log2FoldChange,
    lfc_se_shrunk = shrunken$lfcSE,
    contrast = paste0(numerator, "_vs_", denominator),
    n_numerator = as.integer(sum(model_data$group == numerator)),
    n_denominator = as.integer(sum(model_data$group == denominator)),
    stringsAsFactors = FALSE
  )
}
