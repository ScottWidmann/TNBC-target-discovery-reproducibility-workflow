require_columns <- function(x, required, object_name = deparse(substitute(x))) {
  missing_columns <- setdiff(required, names(x))
  if (length(missing_columns) > 0) {
    stop(
      object_name, " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

discover_expression_files <- function(root) {
  if (!dir.exists(root)) {
    stop("Expression root does not exist: ", root, call. = FALSE)
  }
  paths <- list.files(
    root,
    pattern = "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$|\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(paths) == 0) {
    stop("No STAR count TSV files found under: ", root, call. = FALSE)
  }
  paths <- sort(normalizePath(paths, mustWork = TRUE))
  file_ids <- basename(dirname(paths))
  counts_per_id <- table(file_ids)
  invalid_ids <- names(counts_per_id[counts_per_id != 1L])
  if (length(invalid_ids) > 0) {
    stop(
      "Expected exactly one count file per file UUID directory; invalid: ",
      paste(invalid_ids, collapse = ", "),
      call. = FALSE
    )
  }
  data.frame(
    file_id = file_ids,
    file_name = basename(paths),
    path = paths,
    stringsAsFactors = FALSE
  )
}

validate_expression_metadata <- function(local_files, metadata) {
  require_columns(local_files, "file_id", "local expression files")
  require_columns(metadata, "file_id", "GDC expression metadata")
  duplicate_ids <- unique(metadata$file_id[duplicated(metadata$file_id)])
  if (length(duplicate_ids) > 0) {
    stop(
      "GDC expression metadata contains duplicate file IDs: ",
      paste(duplicate_ids, collapse = ", "),
      call. = FALSE
    )
  }
  missing_ids <- setdiff(local_files$file_id, metadata$file_id)
  if (length(missing_ids) > 0) {
    stop(
      "GDC expression metadata is missing local file IDs: ",
      paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }
  unexpected_ids <- setdiff(metadata$file_id, local_files$file_id)
  if (length(unexpected_ids) > 0) {
    stop(
      "GDC expression metadata contains nonlocal file IDs: ",
      paste(unexpected_ids, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

write_checksum_manifest <- function(paths, output_file = NULL) {
  paths <- sort(normalizePath(paths, mustWork = TRUE))
  checksums <- vapply(
    paths,
    digest::digest,
    algo = "sha256",
    file = TRUE,
    FUN.VALUE = character(1)
  )
  manifest <- data.frame(
    sha256 = unname(checksums),
    path = paths,
    stringsAsFactors = FALSE
  )
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    writeLines(paste(manifest$sha256, manifest$path), output_file, useBytes = TRUE)
  }
  manifest
}

first_existing_column <- function(x, candidates, object_name) {
  selected <- candidates[candidates %in% names(x)]
  if (length(selected) == 0) {
    stop(
      object_name, " lacks all recognized columns: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  selected[[1]]
}

fetch_gdc_file_metadata <- function(local_files, query_fn) {
  require_columns(local_files, c("file_id", "file_name"), "local expression files")
  query_results <- as.data.frame(query_fn(), stringsAsFactors = FALSE)
  id_column <- first_existing_column(query_results, c("file_id", "id"), "GDC query result")
  file_column <- first_existing_column(query_results, c("file_name", "filename"), "GDC query result")
  sample_column <- first_existing_column(
    query_results,
    c("sample_submitter_id", "sample_id", "cases.samples.submitter_id", "cases"),
    "GDC query result"
  )
  type_column <- first_existing_column(
    query_results,
    c("sample_type", "cases.samples.sample_type"),
    "GDC query result"
  )
  patient_candidates <- c("patient_id", "case_submitter_id", "cases.submitter_id")
  patient_column <- patient_candidates[patient_candidates %in% names(query_results)]

  metadata <- data.frame(
    file_id = as.character(query_results[[id_column]]),
    file_name = as.character(query_results[[file_column]]),
    sample_id = as.character(query_results[[sample_column]]),
    sample_type = as.character(query_results[[type_column]]),
    stringsAsFactors = FALSE
  )
  if (length(patient_column) > 0) {
    metadata$patient_id <- as.character(query_results[[patient_column[[1]]]])
  } else {
    metadata$patient_id <- substr(metadata$sample_id, 1, 12)
  }
  metadata <- metadata[metadata$file_id %in% local_files$file_id, , drop = FALSE]
  metadata <- metadata[match(local_files$file_id, metadata$file_id), , drop = FALSE]
  metadata <- metadata[, c("file_id", "file_name", "sample_id", "patient_id", "sample_type")]
  rownames(metadata) <- NULL
  validate_expression_metadata(local_files, metadata)
  metadata
}

load_or_initialize_gdc_metadata <- function(
  local_files,
  metadata_path,
  initialize = FALSE,
  query_fn,
  now_fn = Sys.time
) {
  if (initialize) {
    if (file.exists(metadata_path)) {
      stop(
        "Frozen GDC expression metadata already exists: ", metadata_path,
        call. = FALSE
      )
    }
    metadata <- fetch_gdc_file_metadata(local_files, query_fn)
    metadata$retrieved_at_utc <- format(
      now_fn(),
      "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    )
    dir.create(dirname(metadata_path), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(metadata, metadata_path, sep = "\t")
    return(metadata)
  }

  if (!file.exists(metadata_path)) {
    stop(
      "Frozen GDC expression metadata is absent. Run 00_verify_inputs.R ",
      "with --initialize-gdc-metadata once.",
      call. = FALSE
    )
  }
  metadata <- data.table::fread(metadata_path, data.table = FALSE)
  require_columns(
    metadata,
    c(
      "file_id", "file_name", "sample_id", "patient_id", "sample_type",
      "retrieved_at_utc"
    ),
    "frozen GDC expression metadata"
  )
  validate_expression_metadata(local_files, metadata)
  metadata <- metadata[match(local_files$file_id, metadata$file_id), , drop = FALSE]
  rownames(metadata) <- NULL
  metadata
}

read_analysis_config <- function(path) {
  if (!file.exists(path)) {
    stop("Analysis configuration does not exist: ", path, call. = FALSE)
  }
  yaml::read_yaml(path)
}

query_gdc_expression_metadata <- function(config) {
  query <- TCGAbiolinks::GDCquery(
    project = config$project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )
  TCGAbiolinks::getResults(query)
}
