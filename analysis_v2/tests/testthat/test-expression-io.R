write_star_fixture <- function(path, rows, column_name = "unstranded") {
  header <- paste(c("gene_id", "gene_name", "gene_type", column_name), collapse = "\t")
  body <- apply(rows, 1, paste, collapse = "\t")
  writeLines(c("# gene-model: GENCODE v36", header, body), path)
}

test_that("STAR parser removes summary rows and preserves duplicate symbols", {
  path <- tempfile(fileext = ".tsv")
  rows <- rbind(
    c("N_unmapped", "", "", "10"),
    c("ENSG000001.2", "DUP", "protein_coding", "11"),
    c("ENSG000002.7", "DUP", "lncRNA", "13")
  )
  write_star_fixture(path, rows)

  x <- read_star_count_file(path)

  expect_identical(x$ensembl_id, c("ENSG000001", "ENSG000002"))
  expect_identical(x$ensembl_id_versioned, c("ENSG000001.2", "ENSG000002.7"))
  expect_identical(x$gene_symbol, c("DUP", "DUP"))
  expect_identical(x$count, c(11L, 13L))
})

test_that("STAR parser rejects missing unstranded counts", {
  path <- tempfile(fileext = ".tsv")
  rows <- rbind(c("ENSG000001.2", "A", "protein_coding", "11"))
  write_star_fixture(path, rows, column_name = "tpm_unstranded")

  expect_error(read_star_count_file(path), "unstranded")
})

test_that("STAR parser rejects negative and noninteger counts", {
  negative <- tempfile(fileext = ".tsv")
  fractional <- tempfile(fileext = ".tsv")
  write_star_fixture(negative, rbind(c("ENSG000001.2", "A", "protein_coding", "-1")))
  write_star_fixture(fractional, rbind(c("ENSG000001.2", "A", "protein_coding", "1.5")))

  expect_error(read_star_count_file(negative), "nonnegative integers")
  expect_error(read_star_count_file(fractional), "nonnegative integers")
})

test_that("STAR parser rejects duplicate version-stripped Ensembl IDs", {
  path <- tempfile(fileext = ".tsv")
  rows <- rbind(
    c("ENSG000001.1", "A", "protein_coding", "11"),
    c("ENSG000001.2", "A", "protein_coding", "12")
  )
  write_star_fixture(path, rows)

  expect_error(read_star_count_file(path), "duplicate.*ENSG000001")
})

test_that("matrix assembly uses frozen sample IDs and checks annotations", {
  root <- tempfile("matrix-")
  dir.create(root)
  first <- file.path(root, "first.tsv")
  second <- file.path(root, "second.tsv")
  rows_first <- rbind(
    c("ENSG000001.1", "A", "protein_coding", "11"),
    c("ENSG000002.1", "B", "lncRNA", "12")
  )
  rows_second <- rbind(
    c("ENSG000001.1", "A", "protein_coding", "21"),
    c("ENSG000002.1", "B", "lncRNA", "22")
  )
  write_star_fixture(first, rows_first)
  write_star_fixture(second, rows_second)
  metadata <- data.frame(
    file_id = c("f1", "f2"),
    path = c(first, second),
    sample_id = c("TCGA-AA-0001-01A", "TCGA-AA-0002-11A")
  )

  assembled <- assemble_count_matrix(metadata)

  expect_identical(rownames(assembled$counts), c("ENSG000001", "ENSG000002"))
  expect_identical(colnames(assembled$counts), metadata$sample_id)
  expect_identical(unname(assembled$counts[, 2]), c(21L, 22L))
  expect_identical(assembled$gene_annotation$gene_symbol, c("A", "B"))
})

test_that("matrix assembly rejects changed gene order or annotation", {
  root <- tempfile("matrix-")
  dir.create(root)
  first <- file.path(root, "first.tsv")
  second <- file.path(root, "second.tsv")
  write_star_fixture(
    first,
    rbind(
      c("ENSG000001.1", "A", "protein_coding", "11"),
      c("ENSG000002.1", "B", "lncRNA", "12")
    )
  )
  write_star_fixture(
    second,
    rbind(
      c("ENSG000002.1", "B", "lncRNA", "22"),
      c("ENSG000001.1", "CHANGED", "protein_coding", "21")
    )
  )
  metadata <- data.frame(
    file_id = c("f1", "f2"),
    path = c(first, second),
    sample_id = c("S1", "S2")
  )

  expect_error(assemble_count_matrix(metadata), "gene order or annotation")
})
