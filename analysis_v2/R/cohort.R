normalize_receptor_status <- function(x) {
  value <- trimws(toupper(as.character(x)))
  missing_values <- c(
    "", "NA", "N/A", "[NOT AVAILABLE]", "[NOT APPLICABLE]", "[UNKNOWN]",
    "NOT REPORTED", "UNKNOWN"
  )
  normalized <- rep("Ambiguous", length(value))
  normalized[is.na(x) | is.na(value) | value %in% missing_values] <- "Missing"
  normalized[value == "NEGATIVE"] <- "Negative"
  normalized[value == "POSITIVE"] <- "Positive"
  factor(
    normalized,
    levels = c("Negative", "Positive", "Ambiguous", "Missing")
  )
}

classify_tumor <- function(er, pr, her2) {
  er <- as.character(er)
  pr <- as.character(pr)
  her2 <- as.character(her2)
  all_negative <- er == "Negative" & pr == "Negative" & her2 == "Negative"
  any_positive <- er == "Positive" | pr == "Positive" | her2 == "Positive"
  result <- rep("Unknown_tumor", length(er))
  result[all_negative] <- "TNBC"
  result[!all_negative & any_positive] <- "Other_receptor_defined"
  result
}

read_bcr_patient <- function(path) {
  if (!file.exists(path)) {
    stop("BCR clinical patient file does not exist: ", path, call. = FALSE)
  }
  clinical <- data.table::fread(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    na.strings = c("", "NA", "[Not Available]", "[Not Applicable]", "[Unknown]"),
    data.table = FALSE
  )
  patient_column <- first_existing_column(
    clinical,
    c("bcr_patient_barcode", "submitter_id", "patient_id"),
    "BCR clinical patient table"
  )
  clinical <- clinical[
    grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", clinical[[patient_column]]),
    ,
    drop = FALSE
  ]
  if (nrow(clinical) == 0) {
    stop("BCR clinical patient table contains no TCGA patient rows", call. = FALSE)
  }

  extract_column <- function(candidates, default = NA_character_) {
    present <- candidates[candidates %in% names(clinical)]
    if (length(present) == 0) {
      return(rep(default, nrow(clinical)))
    }
    clinical[[present[[1]]]]
  }
  age <- suppressWarnings(as.numeric(extract_column(
    c("age_at_diagnosis", "age_at_initial_pathologic_diagnosis", "age")
  )))
  result <- data.frame(
    patient_id = as.character(clinical[[patient_column]]),
    er_status = as.character(extract_column(c(
      "er_status_by_ihc", "er_status", "estrogen_receptor_status"
    ))),
    pr_status = as.character(extract_column(c(
      "pr_status_by_ihc", "pr_status", "progesterone_receptor_status"
    ))),
    her2_status = as.character(extract_column(c(
      "her2_status_by_ihc", "her2_status",
      "her2_neu_immunohistochemistry_receptor_status"
    ))),
    age_at_diagnosis = age,
    race = as.character(extract_column(c("race", "race_list", "race_category"))),
    ethnicity = as.character(extract_column(c(
      "ethnicity", "hispanic_or_latino_ethnicity"
    ))),
    tissue_source_site = as.character(extract_column("tissue_source_site")),
    histological_type = as.character(extract_column(c(
      "histological_type", "histologic_diagnosis", "primary_diagnosis"
    ))),
    stringsAsFactors = FALSE
  )
  duplicate_patients <- unique(result$patient_id[duplicated(result$patient_id)])
  if (length(duplicate_patients) > 0) {
    stop(
      "BCR clinical patient table contains duplicate patient IDs: ",
      paste(duplicate_patients, collapse = ", "),
      call. = FALSE
    )
  }
  result
}

build_sample_manifest <- function(metadata, clinical) {
  require_columns(
    metadata,
    c("file_id", "file_name", "path", "sample_id", "patient_id", "sample_type"),
    "expression metadata"
  )
  require_columns(
    clinical,
    c("patient_id", "er_status", "pr_status", "her2_status"),
    "clinical data"
  )
  if (anyDuplicated(metadata$sample_id)) {
    stop("Expression metadata contains duplicate sample IDs", call. = FALSE)
  }
  if (anyDuplicated(clinical$patient_id)) {
    stop("Clinical data contains duplicate patient IDs", call. = FALSE)
  }

  clinical_index <- match(metadata$patient_id, clinical$patient_id)
  result <- metadata
  clinical_columns <- setdiff(names(clinical), "patient_id")
  for (column in clinical_columns) {
    result[[column]] <- clinical[[column]][clinical_index]
  }
  result$er_status <- normalize_receptor_status(result$er_status)
  result$pr_status <- normalize_receptor_status(result$pr_status)
  result$her2_status <- normalize_receptor_status(result$her2_status)

  result$analysis_group <- "Excluded_sample_type"
  is_normal <- result$sample_type == "Solid Tissue Normal"
  is_primary <- result$sample_type == "Primary Tumor"
  result$analysis_group[is_normal] <- "Normal"
  result$analysis_group[is_primary] <- classify_tumor(
    result$er_status[is_primary],
    result$pr_status[is_primary],
    result$her2_status[is_primary]
  )
  result$eligible <- result$analysis_group %in% c(
    "Normal", "TNBC", "Other_receptor_defined"
  )
  result$exclusion_reason <- NA_character_
  result$exclusion_reason[result$analysis_group == "Unknown_tumor"] <-
    "ambiguous_or_missing_receptor"
  result$exclusion_reason[result$analysis_group == "Excluded_sample_type"] <-
    "unsupported_sample_type"
  result
}

select_one_per_patient_group <- function(manifest) {
  require_columns(
    manifest,
    c(
      "patient_id", "sample_id", "analysis_group", "eligible",
      "exclusion_reason"
    ),
    "sample manifest"
  )
  result <- manifest
  result$deduplication_rank <- NA_integer_
  eligible_indices <- which(!is.na(result$eligible) & result$eligible)
  if (length(eligible_indices) > 0) {
    keys <- paste(
      result$patient_id[eligible_indices],
      result$analysis_group[eligible_indices],
      sep = "\r"
    )
    groups <- split(eligible_indices, keys)
    for (indices in groups) {
      ordered <- indices[order(result$sample_id[indices])]
      result$deduplication_rank[ordered] <- seq_along(ordered)
    }
  }
  result$selected_for_analysis <- result$eligible & result$deduplication_rank == 1L
  result$selected_for_analysis[is.na(result$selected_for_analysis)] <- FALSE
  duplicate_rows <- result$eligible & !result$selected_for_analysis
  duplicate_rows[is.na(duplicate_rows)] <- FALSE
  result$exclusion_reason[duplicate_rows] <- "duplicate_patient_group"

  selected <- result[result$selected_for_analysis, , drop = FALSE]
  tnbc_patients <- selected$patient_id[selected$analysis_group == "TNBC"]
  normal_patients <- selected$patient_id[selected$analysis_group == "Normal"]
  paired_patients <- intersect(tnbc_patients, normal_patients)
  result$paired_tnbc_normal <- result$selected_for_analysis &
    result$patient_id %in% paired_patients &
    result$analysis_group %in% c("TNBC", "Normal")
  result
}
