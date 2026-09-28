test_that("ORA background contains only tested and mapped genes", {
  de <- data.frame(
    ensembl_id = c("E1", "E2", "E3"),
    padj = c(0.01, 0.20, NA),
    log2fc_shrunk = c(2, 0, 1)
  )
  terms <- data.frame(
    term_id = c("P", "P", "P"),
    term_name = c("Pathway", "Pathway", "Pathway"),
    ensembl_id = c("E1", "E2", "E9")
  )

  inputs <- prepare_ora_inputs(de, terms, "up", 0.05, 1)

  expect_setequal(inputs$universe, c("E1", "E2"))
  expect_identical(inputs$genes, "E1")
})

test_that("up and down ORA sets use shrunken effect direction", {
  de <- data.frame(
    ensembl_id = c("UP", "DOWN", "SMALL", "NS"),
    padj = c(0.01, 0.01, 0.01, 0.20),
    log2fc_shrunk = c(1.5, -2, 0.5, 3)
  )
  terms <- data.frame(
    term_id = "P",
    term_name = "Pathway",
    ensembl_id = de$ensembl_id
  )

  expect_identical(prepare_ora_inputs(de, terms, "up", 0.05, 1)$genes, "UP")
  expect_identical(prepare_ora_inputs(de, terms, "down", 0.05, 1)$genes, "DOWN")
  expect_error(prepare_ora_inputs(de, terms, "sideways", 0.05, 1), "direction")
})

test_that("mapping audit reports tested, mapped, and unmapped identifiers", {
  de <- data.frame(
    ensembl_id = c("E1", "E2", "E3", "E4"),
    padj = c(0.01, 0.20, NA, 0.5)
  )
  terms <- data.frame(term_id = c("P", "P"), ensembl_id = c("E1", "E4"))

  audit <- map_tested_genes(de, terms, database = "fixture")$audit

  expect_identical(audit$total_de_rows, 4L)
  expect_identical(audit$tested_genes, 3L)
  expect_identical(audit$mapped_tested_genes, 2L)
  expect_identical(audit$unmapped_tested_genes, 1L)
})

test_that("empty ORA foreground returns an explicit status without enrichment rows", {
  de <- data.frame(
    ensembl_id = c("E1", "E2"),
    padj = c(0.5, 0.6),
    log2fc_shrunk = c(0.1, -0.2)
  )
  terms <- data.frame(
    term_id = c("P", "P"),
    term_name = c("Pathway", "Pathway"),
    ensembl_id = c("E1", "E2")
  )

  result <- run_ora(de, terms, "up", 0.05, 1, min_size = 1)

  expect_equal(nrow(result), 0)
  expect_identical(attr(result, "status"), "no_significant_genes")
  expect_identical(attr(result, "foreground_size"), 0L)
  expect_identical(attr(result, "universe_size"), 2L)
})

test_that("ORA results carry foreground and universe sizes", {
  de <- data.frame(
    ensembl_id = paste0("E", 1:10),
    padj = c(rep(0.01, 3), rep(0.5, 7)),
    log2fc_shrunk = c(rep(2, 3), rep(0, 7))
  )
  terms <- data.frame(
    term_id = c(rep("P1", 5), rep("P2", 5)),
    term_name = c(rep("Pathway 1", 5), rep("Pathway 2", 5)),
    ensembl_id = paste0("E", 1:10)
  )

  result <- run_ora(de, terms, "up", 0.05, 1, min_size = 1)

  expect_true(all(c("foreground_size", "universe_size", "direction") %in% names(result)))
  expect_true(all(result$foreground_size == 3L))
  expect_true(all(result$universe_size == 10L))
  expect_true(all(result$direction == "up"))
})

test_that("ranked list uses unique finite statistics sorted decreasing", {
  de <- data.frame(
    ensembl_id = c("E1", "E2", "E3", "E2"),
    stat = c(-2, 3, NA, 1)
  )

  ranked <- prepare_ranked_list(de)

  expect_identical(names(ranked), c("E2", "E1"))
  expect_identical(unname(ranked), c(3, -2))
})
