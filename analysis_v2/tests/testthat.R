library(testthat)

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_argument) != 1) {
  stop("Run tests with Rscript analysis_v2/tests/testthat.R")
}

script_file <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
r_dir <- file.path(analysis_root, "R")

if (dir.exists(r_dir)) {
  r_files <- sort(list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE))
  invisible(lapply(r_files, source))
}

test_dir(file.path(analysis_root, "tests", "testthat"), reporter = "summary")
