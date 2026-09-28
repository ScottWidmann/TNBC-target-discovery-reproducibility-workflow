test_that("expression discovery derives file UUIDs from frozen paths", {
  root <- tempfile("counts-")
  dir.create(file.path(root, "file-uuid"), recursive = TRUE)
  file.create(file.path(root, "file-uuid", "counts.tsv"))

  found <- discover_expression_files(root)

  expect_identical(found$file_id, "file-uuid")
  expect_identical(found$file_name, "counts.tsv")
  expect_identical(found$path, normalizePath(file.path(root, "file-uuid", "counts.tsv")))
})

test_that("expression discovery rejects duplicate files within one UUID directory", {
  root <- tempfile("counts-")
  dir.create(file.path(root, "file-uuid"), recursive = TRUE)
  file.create(file.path(root, "file-uuid", "counts-1.tsv"))
  file.create(file.path(root, "file-uuid", "counts-2.tsv"))

  expect_error(discover_expression_files(root), "exactly one count file")
})

test_that("metadata validation rejects missing and duplicate file IDs", {
  local <- data.frame(
    file_id = c("a", "b"),
    file_name = c("a.tsv", "b.tsv"),
    stringsAsFactors = FALSE
  )

  expect_error(
    validate_expression_metadata(local, data.frame(file_id = "a")),
    "missing.*b"
  )

  duplicate <- data.frame(file_id = c("a", "a", "b"))
  expect_error(validate_expression_metadata(local, duplicate), "duplicate.*a")
})

test_that("required-column errors identify the object and missing fields", {
  x <- data.frame(file_id = "a")

  expect_error(
    require_columns(x, c("file_id", "sample_id", "sample_type"), "expression metadata"),
    "expression metadata.*sample_id.*sample_type"
  )
})

test_that("checksum manifest is deterministic and detects content changes", {
  root <- tempfile("checksums-")
  dir.create(root)
  first <- file.path(root, "a.txt")
  second <- file.path(root, "b.txt")
  writeLines("alpha", first)
  writeLines("beta", second)

  before <- write_checksum_manifest(c(second, first))
  writeLines("changed", first)
  after <- write_checksum_manifest(c(first, second))

  expect_identical(before$path, sort(normalizePath(c(first, second))))
  expect_false(before$sha256[before$path == normalizePath(first)] ==
                 after$sha256[after$path == normalizePath(first)])
  expect_identical(before$sha256[before$path == normalizePath(second)],
                   after$sha256[after$path == normalizePath(second)])
})

test_that("injected metadata query is filtered to local file UUIDs", {
  local <- data.frame(
    file_id = c("a", "b"),
    file_name = c("a.tsv", "b.tsv"),
    stringsAsFactors = FALSE
  )
  query_fn <- function() {
    data.frame(
      id = c("unused", "b", "a"),
      file_name = c("unused.tsv", "b.tsv", "a.tsv"),
      cases = c("TCGA-ZZ-9999", "TCGA-AA-0002", "TCGA-AA-0001"),
      sample_submitter_id = c("TCGA-ZZ-9999-01A", "TCGA-AA-0002-11A", "TCGA-AA-0001-01A"),
      sample_type = c("Primary Tumor", "Solid Tissue Normal", "Primary Tumor"),
      stringsAsFactors = FALSE
    )
  }

  metadata <- fetch_gdc_file_metadata(local, query_fn)

  expect_identical(metadata$file_id, c("a", "b"))
  expect_identical(metadata$patient_id, c("TCGA-AA-0001", "TCGA-AA-0002"))
  expect_identical(metadata$sample_type, c("Primary Tumor", "Solid Tissue Normal"))
})

test_that("offline metadata loading never calls the query function", {
  root <- tempfile("metadata-")
  dir.create(root)
  metadata_path <- file.path(root, "metadata.tsv")
  data.table::fwrite(
    data.frame(
      file_id = "a",
      file_name = "a.tsv",
      sample_id = "TCGA-AA-0001-01A",
      patient_id = "TCGA-AA-0001",
      sample_type = "Primary Tumor",
      retrieved_at_utc = "2026-08-20T00:00:00Z"
    ),
    metadata_path,
    sep = "\t"
  )
  local <- data.frame(file_id = "a", file_name = "a.tsv")
  forbidden_query <- function() stop("network query was called")

  loaded <- load_or_initialize_gdc_metadata(
    local,
    metadata_path,
    initialize = FALSE,
    query_fn = forbidden_query
  )

  expect_identical(loaded$sample_id, "TCGA-AA-0001-01A")
})

test_that("offline metadata loading fails when the frozen table is absent", {
  local <- data.frame(file_id = "a", file_name = "a.tsv")
  missing_path <- tempfile("absent-metadata-")

  expect_error(
    load_or_initialize_gdc_metadata(
      local,
      missing_path,
      initialize = FALSE,
      query_fn = function() stop("must not be called")
    ),
    "initialize-gdc-metadata"
  )
})

test_that("metadata initialization writes retrieval time and refuses overwrite", {
  root <- tempfile("metadata-")
  dir.create(root)
  metadata_path <- file.path(root, "metadata.tsv")
  local <- data.frame(file_id = "a", file_name = "a.tsv")
  query_fn <- function() {
    data.frame(
      id = "a",
      file_name = "a.tsv",
      cases = "TCGA-AA-0001",
      sample_submitter_id = "TCGA-AA-0001-01A",
      sample_type = "Primary Tumor"
    )
  }

  initialized <- load_or_initialize_gdc_metadata(
    local,
    metadata_path,
    initialize = TRUE,
    query_fn = query_fn,
    now_fn = function() as.POSIXct("2026-08-20 12:34:56", tz = "UTC")
  )

  expect_true(file.exists(metadata_path))
  expect_identical(initialized$retrieved_at_utc, "2026-08-20T12:34:56Z")
  expect_error(
    load_or_initialize_gdc_metadata(local, metadata_path, TRUE, query_fn),
    "already exists"
  )
})
