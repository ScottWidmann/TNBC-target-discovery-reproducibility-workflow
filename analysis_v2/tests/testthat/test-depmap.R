test_that("only explicit breast subtype labels enter comparator groups", {
  models <- data.frame(
    ModelID = c("T", "E", "H", "U", "O"),
    CellLineName = c("TN", "ER", "HER2", "Unknown", "Ovary"),
    OncotreeLineage = c("Breast", "Breast", "Breast", "Breast", "Ovary/Fallopian Tube"),
    ModelSubtypeFeatures = c("basal_A TNBC", "luminal ER+", "HER2+", NA, "TNBC")
  )

  groups <- classify_depmap_models(models)

  expect_identical(groups$model_id, c("T", "E", "H"))
  expect_identical(groups$model_group, c("TNBC", "Non_TNBC_breast", "Non_TNBC_breast"))
  expect_identical(groups$source_label, c("basal_A TNBC", "luminal ER+", "HER2+"))
})

test_that("DepMap gene labels remove only trailing Entrez annotations", {
  labels <- c("A1BG (1)", "HLA-DRA (3122)", "GENE.1", "ABC_2", "MELK")

  expect_identical(
    normalize_depmap_gene_symbol(labels),
    c("A1BG", "HLA-DRA", "GENE.1", "ABC_2", "MELK")
  )
})

test_that("effect extraction returns long values for measured panel models", {
  effects <- data.frame(
    ModelID = c("T", "E", "UNUSED"),
    `MELK (9833)` = c(-1, -0.2, 0),
    `HLA-DRA (3122)` = c(-0.5, 0.1, 0),
    check.names = FALSE
  )
  manifest <- data.frame(
    model_id = c("T", "E"),
    model_name = c("TN", "ER"),
    model_group = c("TNBC", "Non_TNBC_breast")
  )

  long <- extract_depmap_effects(effects, manifest)

  expect_equal(nrow(long), 4)
  expect_setequal(long$gene_symbol, c("MELK", "HLA-DRA"))
  expect_equal(long$effect[long$model_id == "T" & long$gene_symbol == "MELK"], -1)
  expect_false("UNUSED" %in% long$model_id)
})

test_that("dependency summaries retain panel size and selectivity direction", {
  effects <- data.frame(
    gene_symbol = rep("G", 4),
    model_name = c("TN1", "MDA-MB-231", "ER1", "ER2"),
    model_group = rep(c("TNBC", "Non_TNBC_breast"), each = 2),
    effect = c(-1, -0.8, -0.2, 0)
  )

  out <- summarize_depmap_by_gene(effects, -0.5)

  expect_equal(out$tnbc_n, 2L)
  expect_equal(out$tnbc_median_effect, -0.9)
  expect_equal(out$tnbc_fraction_below_threshold, 1)
  expect_equal(out$non_tnbc_median_effect, -0.1)
  expect_equal(out$tnbc_minus_non_tnbc, -0.8)
  expect_equal(out$mda_mb_231_effect, -0.8)
})

test_that("common-essential annotation distinguishes false from unavailable", {
  summary <- data.frame(gene_symbol = c("MELK", "OTHER"))
  common <- data.frame(Essentials = c("MELK (9833)"))

  annotated <- annotate_common_essentials(summary, common)

  expect_identical(annotated$common_essential, c(TRUE, FALSE))
})

test_that("single-column common-essential files do not split on spaces", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("Essentials", "AAMP (14)", "MELK (9833)"), path)

  common <- read_common_essential_table(path)

  expect_identical(names(common), "Essentials")
  expect_identical(common$Essentials, c("AAMP (14)", "MELK (9833)"))
})

test_that("DepMap join preserves unique Ensembl keys and duplicate gene symbols", {
  de <- data.frame(
    ensembl_id = c("E1", "E2"),
    gene_symbol = c("DUP", "DUP"),
    padj = c(0.01, 0.02)
  )
  dependency <- data.frame(
    gene_symbol = "DUP",
    tnbc_median_effect = -0.8
  )

  joined <- join_depmap_to_de(de, dependency)

  expect_identical(joined$ensembl_id, c("E1", "E2"))
  expect_identical(joined$tnbc_median_effect, c(-0.8, -0.8))
})
