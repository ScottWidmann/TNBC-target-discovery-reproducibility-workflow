strip_ensembl_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

read_star_count_file <- function(path) {
  if (!file.exists(path)) {
    stop("STAR count file does not exist: ", path, call. = FALSE)
  }
  counts <- tryCatch(
    data.table::fread(path, skip = "gene_id", data.table = FALSE),
    error = function(error) {
      stop("Could not parse STAR count file ", path, ": ", error$message, call. = FALSE)
    }
  )
  require_columns(
    counts,
    c("gene_id", "gene_name", "gene_type", "unstranded"),
    paste0("STAR count file ", path)
  )
  counts <- counts[grepl("^ENSG", counts$gene_id), , drop = FALSE]
  if (nrow(counts) == 0) {
    stop("STAR count file contains no Ensembl gene rows: ", path, call. = FALSE)
  }
  values <- suppressWarnings(as.numeric(counts$unstranded))
  invalid_counts <- !is.finite(values) | values < 0 | values != round(values)
  if (any(invalid_counts)) {
    stop("STAR counts must be nonnegative integers: ", path, call. = FALSE)
  }
  ensembl_ids <- strip_ensembl_version(counts$gene_id)
  duplicate_ids <- unique(ensembl_ids[duplicated(ensembl_ids)])
  if (length(duplicate_ids) > 0) {
    stop(
      "STAR count file has duplicate version-stripped Ensembl IDs: ",
      paste(duplicate_ids, collapse = ", "), " in ", path,
      call. = FALSE
    )
  }
  data.frame(
    ensembl_id = ensembl_ids,
    ensembl_id_versioned = as.character(counts$gene_id),
    gene_symbol = as.character(counts$gene_name),
    gene_type = as.character(counts$gene_type),
    count = as.integer(values),
    stringsAsFactors = FALSE
  )
}

assemble_count_matrix <- function(metadata) {
  require_columns(metadata, c("file_id", "path", "sample_id"), "expression metadata")
  if (anyDuplicated(metadata$file_id)) {
    stop("Expression metadata file_id values must be unique", call. = FALSE)
  }
  if (anyDuplicated(metadata$sample_id)) {
    stop("Expression metadata sample_id values must be unique", call. = FALSE)
  }
  if (nrow(metadata) == 0) {
    stop("Expression metadata contains no files", call. = FALSE)
  }

  first <- read_star_count_file(metadata$path[[1]])
  annotation_columns <- c(
    "ensembl_id", "ensembl_id_versioned", "gene_symbol", "gene_type"
  )
  gene_annotation <- first[, annotation_columns, drop = FALSE]
  count_matrix <- matrix(
    0L,
    nrow = nrow(first),
    ncol = nrow(metadata),
    dimnames = list(first$ensembl_id, as.character(metadata$sample_id))
  )
  count_matrix[, 1] <- first$count

  if (nrow(metadata) > 1) {
    for (index in seq.int(2L, nrow(metadata))) {
      current <- read_star_count_file(metadata$path[[index]])
      current_annotation <- current[, annotation_columns, drop = FALSE]
      if (!identical(gene_annotation, current_annotation)) {
        stop(
          "STAR gene order or annotation differs in file: ",
          metadata$path[[index]],
          call. = FALSE
        )
      }
      count_matrix[, index] <- current$count
      if (index %% 100L == 0L) {
        message("Loaded ", index, " of ", nrow(metadata), " STAR count files")
      }
    }
  }
  storage.mode(count_matrix) <- "integer"
  list(counts = count_matrix, gene_annotation = gene_annotation)
}
