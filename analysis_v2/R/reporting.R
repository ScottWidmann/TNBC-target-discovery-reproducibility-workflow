write_output_manifest <- function(paths, output_file) {
  paths <- sort(unique(normalizePath(paths, mustWork = FALSE)))
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Output manifest contains missing files: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  sizes <- file.info(paths)$size
  empty <- paths[is.na(sizes) | sizes == 0]
  if (length(empty)) {
    stop("Output manifest contains empty files: ", paste(empty, collapse = ", "), call. = FALSE)
  }

  manifest <- data.frame(
    sha256 = vapply(
      paths,
      function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE),
      character(1)
    ),
    path = paths,
    bytes = as.numeric(sizes),
    stringsAsFactors = FALSE
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    paste(manifest$sha256, manifest$path),
    output_file,
    useBytes = TRUE
  )
  manifest
}

validate_final_outputs <- function(root) {
  required <- c(
    "results/candidate_evidence.csv",
    "results/corrected_analysis_report.html",
    "plots/pca_primary.png",
    "plots/volcano_primary.png",
    "plots/concordance_paired_primary.png",
    "plots/concordance_specificity_primary.png",
    "plots/depmap_tnbc_vs_non_tnbc.png",
    "provenance/output_manifest.sha256",
    "provenance/warnings.csv",
    "provenance/session_info.txt"
  )
  paths <- file.path(root, required)
  status <- data.frame(
    relative_path = required,
    exists = file.exists(paths),
    bytes = ifelse(file.exists(paths), file.info(paths)$size, NA_real_),
    stringsAsFactors = FALSE
  )
  invalid <- !status$exists | is.na(status$bytes) | status$bytes == 0
  if (any(invalid)) {
    stop(
      "Required final outputs are missing or empty: ",
      paste(status$relative_path[invalid], collapse = ", "),
      call. = FALSE
    )
  }
  status
}

render_corrected_report <- function(output_file, template = NULL) {
  output_file <- normalizePath(output_file, mustWork = FALSE)
  analysis_root <- dirname(dirname(output_file))
  if (is.null(template)) {
    template <- file.path(analysis_root, "corrected_analysis_report.Rmd")
  }
  if (!file.exists(template)) {
    stop("Report template is missing: ", template, call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  if (rmarkdown::pandoc_available()) {
    rendered <- rmarkdown::render(
      input = template,
      output_file = basename(output_file),
      output_dir = dirname(output_file),
      params = list(analysis_root = analysis_root),
      envir = new.env(parent = globalenv()),
      quiet = TRUE
    )
  } else {
    if (!requireNamespace("markdown", quietly = TRUE)) {
      stop("Neither Pandoc nor the R markdown fallback renderer is available", call. = FALSE)
    }
    temporary_markdown <- tempfile(fileext = ".md")
    on.exit(unlink(temporary_markdown), add = TRUE)
    previous_root <- Sys.getenv("TNBC_ANALYSIS_ROOT", unset = NA_character_)
    Sys.setenv(TNBC_ANALYSIS_ROOT = analysis_root)
    on.exit({
      if (is.na(previous_root)) {
        Sys.unsetenv("TNBC_ANALYSIS_ROOT")
      } else {
        Sys.setenv(TNBC_ANALYSIS_ROOT = previous_root)
      }
    }, add = TRUE)
    render_environment <- new.env(parent = globalenv())
    render_environment$params <- list(analysis_root = analysis_root)
    render_environment$report_analysis_root <- analysis_root
    knitr::knit(
      input = template,
      output = temporary_markdown,
      envir = render_environment,
      quiet = TRUE
    )
    markdown::markdownToHTML(
      file = temporary_markdown,
      output = output_file,
      options = c(
        "+toc", "+embed_resources", "+table", "+auto_identifiers",
        "+autolink"
      ),
      title = "Corrected TCGA-BRCA TNBC Target-Nomination Analysis"
    )
    rendered <- output_file
  }
  normalizePath(rendered, mustWork = TRUE)
}

validate_rendered_report <- function(path, minimum_images = 1L) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Rendered report is missing or empty: ", path, call. = FALSE)
  }
  html <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (grepl("## Error", html, fixed = TRUE)) {
    stop("Rendered report contains an execution error", call. = FALSE)
  }
  required_sections <- c(
    "Scope and interpretation",
    "Candidate evidence classes",
    "Warnings and limitations"
  )
  missing_sections <- required_sections[!vapply(
    required_sections,
    grepl,
    logical(1),
    x = html,
    fixed = TRUE
  )]
  if (length(missing_sections)) {
    stop(
      "Rendered report lacks required sections: ",
      paste(missing_sections, collapse = ", "),
      call. = FALSE
    )
  }
  image_matches <- gregexpr("data:image/", html, fixed = TRUE)[[1]]
  embedded_images <- if (identical(image_matches[[1]], -1L)) 0L else length(image_matches)
  if (embedded_images < minimum_images) {
    stop(
      "Rendered report contains ", embedded_images,
      " embedded images; expected at least ", minimum_images,
      call. = FALSE
    )
  }
  data.frame(
    path = normalizePath(path, mustWork = TRUE),
    bytes = file.info(path)$size,
    embedded_images = embedded_images,
    stringsAsFactors = FALSE
  )
}

compute_pca_scores <- function(counts, manifest, n_top = 500) {
  if (is.null(colnames(counts))) {
    stop("Count matrix must have sample column names", call. = FALSE)
  }
  require_columns(manifest, c("sample_id", "analysis_group"), "sample manifest")
  if (!all(colnames(counts) %in% manifest$sample_id)) {
    stop("Count matrix samples are missing from the sample manifest", call. = FALSE)
  }
  library_size <- colSums(counts)
  if (any(!is.finite(library_size)) || any(library_size <= 0)) {
    stop("Count matrix contains a sample with invalid library size", call. = FALSE)
  }
  log_cpm <- log2(sweep(counts, 2, library_size / 1e6, "/") + 1)
  variances <- apply(log_cpm, 1, stats::var)
  keep <- head(order(variances, decreasing = TRUE), min(n_top, length(variances)))
  fit <- stats::prcomp(t(log_cpm[keep, , drop = FALSE]), center = TRUE, scale. = FALSE)
  if (ncol(fit$x) < 2) {
    stop("PCA produced fewer than two components", call. = FALSE)
  }
  scores <- data.frame(
    sample_id = rownames(fit$x),
    PC1 = fit$x[, 1],
    PC2 = fit$x[, 2],
    stringsAsFactors = FALSE
  )
  index <- match(scores$sample_id, manifest$sample_id)
  scores$analysis_group <- manifest$analysis_group[index]
  variance_percent <- 100 * fit$sdev[1:2]^2 / sum(fit$sdev^2)
  list(scores = scores, variance_percent = variance_percent)
}

wrap_plot_caption <- function(x, width = 82) {
  paste(strwrap(x, width = width), collapse = "\n")
}

make_pca_plot <- function(pca, group_counts) {
  caption <- wrap_plot_caption(paste0(
    "Primary cohort, TNBC and Normal; library-size-normalized log2 counts; ",
    "top variable genes; unsupervised display with no significance threshold. ",
    paste(names(group_counts), "n=", as.integer(group_counts), collapse = "; ")
  ))
  ggplot2::ggplot(
    pca$scores,
    ggplot2::aes(x = PC1, y = PC2, color = analysis_group)
  ) +
    ggplot2::geom_point(alpha = 0.7, size = 1.8) +
    ggplot2::labs(
      title = "Global expression structure",
      x = sprintf("PC1 (%.1f%%)", pca$variance_percent[[1]]),
      y = sprintf("PC2 (%.1f%%)", pca$variance_percent[[2]]),
      color = "Analysis group",
      caption = caption
    ) +
    ggplot2::theme_bw(base_size = 11)
}

make_volcano_plot <- function(de, alpha, lfc_threshold) {
  require_columns(
    de,
    c("log2fc_shrunk", "padj", "gene_symbol", "n_numerator", "n_denominator"),
    "primary DE results"
  )
  status <- rep("Not significant", nrow(de))
  status[!is.na(de$padj) & de$padj < alpha & de$log2fc_shrunk > lfc_threshold] <- "Up"
  status[!is.na(de$padj) & de$padj < alpha & de$log2fc_shrunk < -lfc_threshold] <- "Down"
  plot_data <- data.frame(
    log2fc = de$log2fc_shrunk,
    minus_log10_fdr = -log10(pmax(de$padj, .Machine$double.xmin)),
    status = factor(status, levels = c("Down", "Not significant", "Up")),
    stringsAsFactors = FALSE
  )
  n_numerator <- unique(de$n_numerator[!is.na(de$n_numerator)])
  n_denominator <- unique(de$n_denominator[!is.na(de$n_denominator)])
  caption <- wrap_plot_caption(sprintf(
    "Contrast: TNBC - Normal (TNBC n=%s; Normal n=%s). FDR < %.2f and |shrunken log2FC| > %.1f.",
    paste(n_numerator, collapse = "/"),
    paste(n_denominator, collapse = "/"),
    alpha,
    lfc_threshold
  ))
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = log2fc, y = minus_log10_fdr, color = status)
  ) +
    ggplot2::geom_point(alpha = 0.45, size = 0.8, na.rm = TRUE) +
    ggplot2::geom_vline(
      xintercept = c(-lfc_threshold, lfc_threshold),
      linetype = "dashed",
      color = "grey45"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      color = "grey45"
    ) +
    ggplot2::scale_color_manual(values = c(Down = "#2C7BB6", `Not significant` = "grey70", Up = "#D7191C")) +
    ggplot2::labs(
      title = "Primary differential expression",
      x = "Shrunken log2 fold change (TNBC - Normal)",
      y = "-log10 adjusted p-value",
      color = NULL,
      caption = caption
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "top")
}

make_concordance_plot <- function(primary, sensitivity, sensitivity_name) {
  require_columns(primary, c("ensembl_id", "log2fc_shrunk"), "primary DE results")
  require_columns(
    sensitivity,
    c("ensembl_id", "log2fc_shrunk", "n_numerator", "n_denominator"),
    paste(sensitivity_name, "DE results")
  )
  index <- match(primary$ensembl_id, sensitivity$ensembl_id)
  plot_data <- data.frame(
    primary = primary$log2fc_shrunk,
    sensitivity = sensitivity$log2fc_shrunk[index]
  )
  plot_data <- plot_data[stats::complete.cases(plot_data), , drop = FALSE]
  correlation <- if (nrow(plot_data) > 2) stats::cor(plot_data$primary, plot_data$sensitivity) else NA_real_
  n_num <- unique(sensitivity$n_numerator[!is.na(sensitivity$n_numerator)])
  n_den <- unique(sensitivity$n_denominator[!is.na(sensitivity$n_denominator)])
  direction <- if (identical(sensitivity_name, "Paired")) {
    "paired TNBC - paired Normal"
  } else {
    "TNBC - receptor-defined other tumor"
  }
  caption <- wrap_plot_caption(sprintf(
    "Sensitivity contrast: %s (n=%s vs n=%s); primary direction: TNBC - Normal; %s genes with finite effects; no significance filter; Pearson r=%.3f.",
    direction,
    paste(n_num, collapse = "/"),
    paste(n_den, collapse = "/"),
    format(nrow(plot_data), big.mark = ","),
    correlation
  ))
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = primary, y = sensitivity)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75") +
    ggplot2::geom_vline(xintercept = 0, color = "grey75") +
    ggplot2::geom_point(alpha = 0.25, size = 0.7, color = "#2C7FB8") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey35") +
    ggplot2::labs(
      title = paste(sensitivity_name, "effect-direction concordance"),
      x = "Primary shrunken log2FC",
      y = paste(sensitivity_name, "shrunken log2FC"),
      caption = caption
    ) +
    ggplot2::theme_bw(base_size = 11)
}

make_enrichment_dotplot <- function(
    x,
    title,
    max_terms = 20,
    n_tnbc = NA_integer_,
    n_normal = NA_integer_,
    alpha = 0.05,
    lfc_threshold = 1) {
  if (!nrow(x)) {
    return(NULL)
  }
  require_columns(x, c("term_name", "padj", "gene_count"), paste(title, "enrichment"))
  plot_data <- x[!is.na(x$padj) & x$padj < alpha, , drop = FALSE]
  if (!nrow(plot_data)) {
    return(NULL)
  }
  plot_data <- head(plot_data[order(plot_data$padj, -plot_data$gene_count), , drop = FALSE], max_terms)
  plot_data$term_name <- vapply(
    as.character(plot_data$term_name),
    function(term) paste(strwrap(term, width = 42), collapse = "\n"),
    character(1)
  )
  plot_data$term_name <- factor(plot_data$term_name, levels = rev(plot_data$term_name))
  direction <- if (grepl("downregulated", title, ignore.case = TRUE)) "down" else "up"
  caption <- wrap_plot_caption(sprintf(
    "Primary TNBC - Normal contrast (TNBC n=%s; Normal n=%s); %s foreground: FDR < %.2f and shrunken log2FC %s %.1f; universe is the tested, locally mapped gene set.",
    ifelse(is.na(n_tnbc), "not supplied", n_tnbc),
    ifelse(is.na(n_normal), "not supplied", n_normal),
    direction,
    alpha,
    ifelse(direction == "up", ">", "<"),
    ifelse(direction == "up", lfc_threshold, -lfc_threshold)
  ))
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = -log10(pmax(padj, .Machine$double.xmin)), y = term_name, size = gene_count, color = padj)
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c(direction = -1) +
    ggplot2::labs(
      title = title,
      x = "-log10 adjusted p-value",
      y = NULL,
      size = "Genes",
      color = "FDR",
      caption = caption
    ) +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(
        title.position = "top",
        barwidth = grid::unit(6, "cm")
      ),
      size = ggplot2::guide_legend(title.position = "top")
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = ggplot2::element_text(size = 14),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0.5)
    )
}

make_depmap_plot <- function(depmap, candidate_evidence, config, n_tnbc, n_comparator) {
  require_columns(
    depmap,
    c("gene_symbol", "tnbc_median_effect", "non_tnbc_median_effect"),
    "DepMap summary"
  )
  classes <- candidate_evidence$evidence_class[match(
    toupper(depmap$gene_symbol),
    toupper(candidate_evidence$gene_symbol)
  )]
  highlight <- ifelse(classes %in% c("A_multistream", "B_expression_dependency"), classes, "Other")
  plot_data <- data.frame(
    tnbc = depmap$tnbc_median_effect,
    comparator = depmap$non_tnbc_median_effect,
    evidence_class = factor(highlight, levels = c("Other", "B_expression_dependency", "A_multistream"))
  )
  caption <- wrap_plot_caption(sprintf(
    "Median gene effect across TNBC n=%d and ER+/HER2+ breast comparator n=%d models. Dependency <= %.1f; TNBC-minus-comparator selectivity <= %.1f.",
    n_tnbc,
    n_comparator,
    config$depmap$dependency_threshold,
    config$depmap$selectivity_delta
  ))
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = comparator, y = tnbc, color = evidence_class)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = config$depmap$dependency_threshold, linetype = "dotted") +
    ggplot2::geom_point(alpha = 0.45, size = 0.9, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c(Other = "grey70", B_expression_dependency = "#FDAE61", A_multistream = "#D7191C"),
      labels = c(Other = "Other", B_expression_dependency = "Class B", A_multistream = "Class A")
    ) +
    ggplot2::labs(
      title = "DepMap dependency across breast model panels",
      x = "Median gene effect: ER+/HER2+ breast models",
      y = "Median gene effect: TNBC models",
      color = "Evidence class",
      caption = caption
    ) +
    ggplot2::theme_bw(base_size = 11)
}

save_report_plot <- function(plot, path, width = 8, height = 6) {
  if (is.null(plot)) {
    return(FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 180,
    bg = "white"
  )
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Plot was not written: ", path, call. = FALSE)
  }
  TRUE
}
