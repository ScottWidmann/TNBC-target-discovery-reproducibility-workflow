test_that("stage runner records successful boundaries in order", {
  stage_one <- tempfile(fileext = ".R")
  stage_two <- tempfile(fileext = ".R")
  log_path <- tempfile(fileext = ".csv")
  writeLines("quit(status = 0)", stage_one)
  writeLines("quit(status = 0)", stage_two)

  log <- run_analysis_stages(c(stage_one, stage_two), log_path)

  expect_identical(log$stage, basename(c(stage_one, stage_two)))
  expect_identical(log$status, c(0L, 0L))
  expect_true(all(nzchar(log$started_at_utc)))
  expect_true(file.exists(log_path))
})

test_that("stage runner fails fast and preserves the failure record", {
  failing_stage <- tempfile(fileext = ".R")
  unreachable_stage <- tempfile(fileext = ".R")
  log_path <- tempfile(fileext = ".csv")
  writeLines("quit(status = 7)", failing_stage)
  writeLines("quit(status = 0)", unreachable_stage)

  expect_error(
    run_analysis_stages(c(failing_stage, unreachable_stage), log_path),
    "failed with status 7"
  )
  log <- data.table::fread(log_path, data.table = FALSE)
  expect_equal(nrow(log), 1)
  expect_identical(log$status, 7L)
})
