format_utc_time <- function(x) {
  format(as.POSIXct(x), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

run_analysis_stages <- function(
    stages,
    log_path,
    rscript = Sys.which("Rscript"),
    time_fn = Sys.time) {
  if (!nzchar(rscript)) {
    stop("Rscript executable was not found", call. = FALSE)
  }
  missing <- stages[!file.exists(stages)]
  if (length(missing)) {
    stop("Analysis stages are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  log <- data.frame(
    stage = character(),
    command = character(),
    started_at_utc = character(),
    ended_at_utc = character(),
    status = integer(),
    stringsAsFactors = FALSE
  )

  for (stage in stages) {
    stage <- normalizePath(stage, mustWork = TRUE)
    started <- time_fn()
    message("Starting ", basename(stage))
    status <- suppressWarnings(system2(rscript, args = stage))
    ended <- time_fn()
    record <- data.frame(
      stage = basename(stage),
      command = paste(shQuote(rscript), shQuote(stage)),
      started_at_utc = format_utc_time(started),
      ended_at_utc = format_utc_time(ended),
      status = as.integer(status),
      stringsAsFactors = FALSE
    )
    log <- rbind(log, record)
    data.table::fwrite(log, log_path)
    if (status != 0) {
      stop(
        basename(stage), " failed with status ", status,
        "; see ", log_path,
        call. = FALSE
      )
    }
    message("Completed ", basename(stage))
  }
  log
}
