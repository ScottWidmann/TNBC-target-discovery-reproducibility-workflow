test_that("contrast filtering ignores samples outside the contrast", {
  counts <- matrix(
    c(10, 0, 0, 0, 0, 0, 999, 999),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("ENSG1", "ENSG2"), c("a", "b", "outside1", "outside2"))
  )

  kept <- filter_counts_for_contrast(counts, c("a", "b"), 10, 0.5)

  expect_identical(rownames(kept), "ENSG1")
  expect_identical(colnames(kept), c("a", "b"))
})

test_that("contrast filtering rejects absent samples and invalid thresholds", {
  counts <- matrix(10L, nrow = 1, dimnames = list("ENSG1", "a"))

  expect_error(filter_counts_for_contrast(counts, "missing", 10, 0.5), "missing")
  expect_error(filter_counts_for_contrast(counts, "a", -1, 0.5), "min_count")
  expect_error(filter_counts_for_contrast(counts, "a", 10, 0), "min_fraction")
})

test_that("gene joins retain Ensembl primary keys and duplicate symbols", {
  de <- data.frame(ensembl_id = c("ENSG1", "ENSG2"), padj = c(0.01, 0.02))
  genes <- data.frame(
    ensembl_id = c("ENSG1", "ENSG2"),
    gene_symbol = c("DUP", "DUP"),
    gene_type = c("protein_coding", "lncRNA")
  )

  joined <- join_gene_annotation(de, genes)

  expect_identical(joined$ensembl_id, c("ENSG1", "ENSG2"))
  expect_identical(joined$gene_symbol, c("DUP", "DUP"))
  expect_identical(joined$gene_type, c("protein_coding", "lncRNA"))
})

test_that("covariate gate accepts variation and rejects missingness or confounding", {
  col_data <- data.frame(
    group = rep(c("Normal", "TNBC"), each = 5),
    age = c(40, 45, 50, 55, 60, 42, 47, 52, 57, 62),
    site = rep(c("A", "B"), each = 5),
    incomplete = c(1:5, rep(NA, 5))
  )

  audit <- assess_covariates(
    col_data,
    c("age", "site", "incomplete", "absent"),
    min_complete = 0.9
  )

  expect_true(audit$eligible[audit$covariate == "age"])
  expect_match(audit$reason[audit$covariate == "site"], "confounded")
  expect_match(audit$reason[audit$covariate == "incomplete"], "completeness")
  expect_match(audit$reason[audit$covariate == "absent"], "absent")
})

test_that("DE designs place patient before group only for paired models", {
  expect_identical(
    paste(deparse(build_deseq_design(c("age"), paired = FALSE)), collapse = ""),
    "~age + group"
  )
  expect_identical(
    paste(deparse(build_deseq_design(character(), paired = TRUE)), collapse = ""),
    "~patient_id + group"
  )
  expect_error(build_deseq_design("age", paired = TRUE), "Paired models")
})

test_that("numeric model covariates are centered and scaled", {
  col_data <- data.frame(
    age = c(40, 50, 60, 70),
    race = c("A", "A", "B", "B")
  )

  prepared <- prepare_model_covariates(col_data, c("age", "race"))

  expect_equal(mean(prepared$age), 0, tolerance = 1e-12)
  expect_equal(stats::sd(prepared$age), 1, tolerance = 1e-12)
  expect_true(is.factor(prepared$race))
})

test_that("patient blocks preserve equality with safe factor labels", {
  patient_id <- c("TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0001")

  block <- encode_patient_block(patient_id)

  expect_identical(block[[1]], block[[3]])
  expect_false(block[[1]] == block[[2]])
  expect_true(all(grepl("^[A-Za-z0-9_.]+$", levels(block))))
})

test_that("DESeq2 result contract retains raw and shrunken effects by Ensembl ID", {
  sample_ids <- paste0("S", 1:8)
  set.seed(20260820)
  counts <- matrix(
    stats::rnbinom(100 * 8, mu = 50, size = 5),
    nrow = 100,
    ncol = 8
  )
  rownames(counts) <- paste0("ENSG", seq_len(nrow(counts)))
  colnames(counts) <- sample_ids
  counts[1, 1:4] <- c(10, 12, 9, 11)
  counts[1, 5:8] <- c(100, 110, 95, 105)
  storage.mode(counts) <- "integer"
  col_data <- data.frame(
    sample_id = sample_ids,
    group = rep(c("Normal", "TNBC"), each = 4),
    patient_id = paste0("P", 1:8),
    row.names = sample_ids
  )

  result <- run_deseq_contrast(
    counts,
    col_data,
    numerator = "TNBC",
    denominator = "Normal",
    covariates = character(),
    paired = FALSE,
    shrink_type = "apeglm"
  )

  expect_identical(result$ensembl_id, rownames(counts))
  expect_true(all(c(
    "base_mean", "log2fc_raw", "lfc_se_raw", "stat", "pvalue", "padj",
    "log2fc_shrunk", "lfc_se_shrunk", "contrast", "n_numerator", "n_denominator"
  ) %in% names(result)))
  expect_identical(unique(result$contrast), "TNBC_vs_Normal")
  expect_identical(unique(result$n_numerator), 4L)
  expect_gt(result$log2fc_shrunk[result$ensembl_id == "ENSG1"], 0)
})
