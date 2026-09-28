source(testthat::test_path("..", "..", "DepMapCrispr", "R", "depmap_integration.R"))

test_that("DE gene_symbol values survive standardization and filtering", {
  de <- data.frame(
    base_mean = c(200, 300, 400),
    log2_fold_change = c(2.0, 0.5, 1.5),
    padj = c(0.01, 0.01, 0.20),
    stat = c(5, 2, 1),
    gene_symbol = c("MELK", "TP53", "SLC7A5")
  )

  standardized <- standardize_de_table(de)
  filtered <- filter_upregulated_de(standardized)

  expect_identical(standardized$GENE, c("MELK", "TP53", "SLC7A5"))
  expect_identical(filtered$GENE, "MELK")
})

test_that("make.unique suffixes do not duplicate a DepMap gene", {
  de <- data.frame(
    base_mean = c(500, 300),
    log2_fold_change = c(2.0, 1.5),
    padj = c(0.001, 0.01),
    stat = c(7, 5),
    gene_symbol = c("MELK", "MELK.1")
  )

  filtered <- filter_upregulated_de(standardize_de_table(de))

  expect_identical(filtered$GENE, "MELK")
  expect_equal(filtered$padj, 0.001)
})

test_that("DepMap column labels normalize to HUGO symbols", {
  labels <- c("A1BG (1)", "BRCA1 (672)", "tp53_7157", "MELK")

  expect_identical(
    normalize_gene_symbol(labels),
    c("A1BG", "BRCA1", "TP53", "MELK")
  )
})

test_that("an integration returns matched genes without probability data", {
  de <- data.frame(
    base_mean = c(500, 600),
    log2_fold_change = c(2.5, 2.0),
    padj = c(0.001, 0.002),
    stat = c(7, 6),
    gene_symbol = c("MELK", "SLC7A5")
  )
  effects <- data.frame(
    model_id = "ACH-TEST",
    `MELK (9833)` = -0.8,
    `SLC7A5 (8140)` = -0.2,
    check.names = FALSE
  )

  joined <- join_de_with_model_effects(de, effects, "ACH-TEST")
  plot_data <- prepare_plot_data(joined, dependency_threshold = -0.5)

  expect_identical(joined$GENE, c("MELK", "SLC7A5"))
  expect_equal(joined$DEPENDENCY, c(-0.8, -0.2))
  expect_identical(plot_data$hit_base, c(TRUE, FALSE))
  expect_false("P_DEP" %in% names(plot_data))
})

test_that("DepMap extraction preserves punctuation in HUGO symbols", {
  de <- data.frame(
    base_mean = 250,
    log2_fold_change = 1.8,
    padj = 0.01,
    stat = 5,
    gene_symbol = "HLA-DRA"
  )
  effects <- data.frame(
    model_id = "ACH-TEST",
    `HLA-DRA (3122)` = -0.7,
    check.names = FALSE
  )

  joined <- join_de_with_model_effects(de, effects, "ACH-TEST")

  expect_identical(joined$GENE, "HLA-DRA")
  expect_equal(joined$DEPENDENCY, -0.7)
})

test_that("ScreenGeneDependency is recognized as a probability input", {
  input_dir <- tempfile("depmap-inputs-")
  dir.create(input_dir)
  file.create(file.path(input_dir, "ScreenGeneDependency.csv"))

  found <- find_first_existing(
    input_dir,
    c("CRISPRGeneDependency.csv", "ScreenGeneDependency.csv")
  )

  expect_identical(basename(found), "ScreenGeneDependency.csv")
})

test_that("the integration writes nonempty tables and a plot", {
  input_dir <- tempfile("depmap-run-")
  output_dir <- file.path(input_dir, "outputs")
  dir.create(input_dir)

  utils::write.csv(
    data.frame(
      baseMean = c(500, 600),
      log2FoldChange = c(2.5, 2.0),
      padj = c(0.001, 0.002),
      stat = c(7, 6),
      gene_symbol = c("MELK", "SLC7A5")
    ),
    file.path(input_dir, "TNBC_vs_Normal_DEG.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(
      ModelID = "ACH-TEST",
      CellLineName = "MDA-MB-231",
      StrippedCellLineName = "MDAMB231"
    ),
    file.path(input_dir, "Model.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(
      V1 = "ACH-TEST",
      `MELK (9833)` = -0.8,
      `SLC7A5 (8140)` = -0.2,
      check.names = FALSE
    ),
    file.path(input_dir, "CRISPRGeneEffect.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(Essentials = "MELK (9833)"),
    file.path(input_dir, "CRISPRInferredCommonEssentials.csv"),
    row.names = FALSE
  )

  result <- run_depmap_integration(input_dir, output_dir)

  expect_equal(nrow(result$joined), 2)
  expect_identical(result$candidates$GENE, "MELK")
  expect_true(result$candidates$COMMON_ESSENTIAL)
  expect_equal(nrow(result$noncommon_candidates), 0)
  expect_true(file.exists(file.path(output_dir, "DE_up_vs_DepMap_MDA_MB_231_joined.csv")))
  expect_true(file.exists(file.path(output_dir, "Candidates_overexpressed_and_essential_MDA_MB_231.csv")))
  expect_true(file.exists(file.path(output_dir, "Candidates_noncommon_overexpressed_and_essential_MDA_MB_231.csv")))
  plot_file <- file.path(output_dir, "Scatter_log2FC_vs_Dependency_MDA_MB_231.png")
  expect_true(file.exists(plot_file))

  first_plot_hash <- unname(tools::md5sum(plot_file))
  run_depmap_integration(input_dir, output_dir)
  expect_identical(unname(tools::md5sum(plot_file)), first_plot_hash)
})
