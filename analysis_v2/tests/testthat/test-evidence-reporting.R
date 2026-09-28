test_that("candidate evidence is classified without an arithmetic score", {
  x <- data.frame(
    primary_hit = c(TRUE, TRUE, TRUE, FALSE),
    paired_concordant = c(TRUE, FALSE, FALSE, TRUE),
    specificity_concordant = c(TRUE, TRUE, FALSE, TRUE),
    tnbc_dependency = c(TRUE, TRUE, FALSE, TRUE),
    depmap_selective = c(TRUE, FALSE, FALSE, TRUE),
    common_essential = c(FALSE, FALSE, FALSE, FALSE)
  )

  out <- assign_evidence_class(x, list())

  expect_identical(
    out,
    c("A_multistream", "B_expression_dependency", "C_expression_only", "not_candidate")
  )
  expect_false(any(grepl("score", names(x), ignore.case = TRUE)))
})

test_that("common essentials cannot enter the multistream class", {
  x <- data.frame(
    primary_hit = TRUE,
    paired_concordant = TRUE,
    specificity_concordant = TRUE,
    tnbc_dependency = TRUE,
    depmap_selective = TRUE,
    common_essential = TRUE
  )

  expect_identical(assign_evidence_class(x, list()), "B_expression_dependency")
})

test_that("evidence assembly applies declared thresholds and records unavailable streams", {
  primary <- data.frame(
    ensembl_id = c("E1", "E2", "E3"),
    gene_symbol = c("A", "B", "C"),
    gene_type = "protein_coding",
    padj = c(0.01, 0.01, 0.20),
    log2fc_shrunk = c(1.5, 1.2, 2.0)
  )
  paired <- data.frame(
    ensembl_id = c("E1", "E2"),
    padj = c(0.1, NA),
    log2fc_shrunk = c(0.2, NA)
  )
  specificity <- data.frame(
    ensembl_id = c("E1", "E2", "E3"),
    padj = c(0.01, 0.01, 0.01),
    log2fc_shrunk = c(0.5, -0.4, 0.5)
  )
  depmap <- data.frame(
    gene_symbol = c("A", "B"),
    tnbc_median_effect = c(-0.8, -0.8),
    tnbc_fraction_below_threshold = c(0.8, 0.4),
    tnbc_minus_non_tnbc = c(-0.4, -0.4),
    common_essential = c(FALSE, FALSE)
  )
  tractability <- data.frame(
    ensembl_id = character(),
    gene_symbol = character(),
    tractability_status = character(),
    source_name = character(),
    source_version = character(),
    source_record_id = character(),
    evidence_note = character()
  )
  config <- list(
    differential_expression = list(alpha = 0.05, lfc_threshold = 1),
    depmap = list(
      dependency_threshold = -0.5,
      dependency_fraction = 0.5,
      selectivity_delta = -0.2
    )
  )

  out <- assemble_candidate_evidence(
    primary, paired, specificity, depmap, tractability, config
  )

  expect_identical(out$evidence_class, c("A_multistream", "C_expression_only", "not_candidate"))
  expect_identical(out$paired_available, c(TRUE, FALSE, FALSE))
  expect_identical(out$paired_concordant, c(TRUE, FALSE, FALSE))
  expect_identical(out$tnbc_dependency, c(TRUE, FALSE, FALSE))
  expect_true(all(out$tractability_assessment == "not_assessed"))
  expect_false(any(grepl("score", names(out), ignore.case = TRUE)))
})

test_that("evidence assembly preserves duplicate symbols and unique Ensembl keys", {
  primary <- data.frame(
    ensembl_id = c("E1", "E2"),
    gene_symbol = c("DUP", "DUP"),
    gene_type = c("a", "b"),
    padj = c(0.01, 0.02),
    log2fc_shrunk = c(2, 2)
  )
  sensitivity <- data.frame(
    ensembl_id = c("E1", "E2"),
    padj = c(0.01, 0.01),
    log2fc_shrunk = c(1, 1)
  )
  depmap <- data.frame(
    gene_symbol = "DUP",
    tnbc_median_effect = -1,
    tnbc_fraction_below_threshold = 1,
    tnbc_minus_non_tnbc = -0.5,
    common_essential = FALSE
  )
  tractability <- data.frame(
    ensembl_id = character(), gene_symbol = character(),
    tractability_status = character(), source_name = character(),
    source_version = character(), source_record_id = character(),
    evidence_note = character()
  )
  config <- list(
    differential_expression = list(alpha = 0.05, lfc_threshold = 1),
    depmap = list(dependency_threshold = -0.5, dependency_fraction = 0.5, selectivity_delta = -0.2)
  )

  out <- assemble_candidate_evidence(
    primary, sensitivity, sensitivity, depmap, tractability, config
  )

  expect_identical(out$ensembl_id, c("E1", "E2"))
  expect_identical(out$tnbc_median_effect, c(-1, -1))
})

test_that("output manifest rejects missing or empty claimed files", {
  root <- tempfile("outputs-")
  dir.create(root)
  good <- file.path(root, "good.csv")
  empty <- file.path(root, "empty.csv")
  writeLines("x\n1", good)
  file.create(empty)

  expect_error(
    write_output_manifest(c(good, file.path(root, "missing.csv")), tempfile()),
    "missing"
  )
  expect_error(write_output_manifest(c(good, empty), tempfile()), "empty")
})

test_that("PCA scores retain sample groups and report explained variance", {
  counts <- matrix(
    c(10, 12, 100, 110, 20, 22, 200, 220, 5, 6, 50, 55),
    nrow = 3,
    dimnames = list(paste0("E", 1:3), paste0("S", 1:4))
  )
  manifest <- data.frame(
    sample_id = paste0("S", 1:4),
    analysis_group = c("Normal", "Normal", "TNBC", "TNBC")
  )

  pca <- compute_pca_scores(counts, manifest, n_top = 3)

  expect_identical(pca$scores$sample_id, paste0("S", 1:4))
  expect_identical(pca$scores$analysis_group, manifest$analysis_group)
  expect_length(pca$variance_percent, 2)
  expect_true(all(is.finite(pca$variance_percent)))
})

test_that("primary volcano caption declares direction, sample sizes, and thresholds", {
  de <- data.frame(
    log2fc_shrunk = c(2, -2, 0),
    padj = c(0.01, 0.02, 0.5),
    gene_symbol = c("UP", "DOWN", "NONE"),
    n_numerator = 115,
    n_denominator = 113
  )

  plot <- make_volcano_plot(de, alpha = 0.05, lfc_threshold = 1)

  expect_s3_class(plot, "ggplot")
  expect_match(plot$labels$caption, "TNBC - Normal", fixed = TRUE)
  expect_match(plot$labels$caption, "n=115", fixed = TRUE)
  expect_match(plot$labels$caption, "FDR < 0.05", fixed = TRUE)
})

test_that("empty enrichment results yield an explicit no-plot status", {
  empty <- data.frame(term_name = character(), padj = numeric(), gene_count = integer())

  expect_null(make_enrichment_dotplot(empty, "GO BP up"))
})

test_that("rendered-report validation rejects embedded execution errors", {
  bad <- tempfile(fileext = ".html")
  good <- tempfile(fileext = ".html")
  writeLines("<html><body><pre>## Error in chunk: missing plot</pre></body></html>", bad)
  writeLines(
    paste0(
      "<html><body>",
      "Scope and interpretation Candidate evidence classes ",
      "Warnings and limitations <img src='data:image/png;base64,AA=='>",
      "</body></html>"
    ),
    good
  )

  expect_error(validate_rendered_report(bad), "execution error")
  expect_true(validate_rendered_report(good)$embedded_images >= 1)
})
