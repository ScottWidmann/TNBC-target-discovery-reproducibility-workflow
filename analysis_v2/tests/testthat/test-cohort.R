test_that("receptor normalization distinguishes negative, positive, ambiguous, and missing", {
  observed <- normalize_receptor_status(c(
    "Negative", "positive", "Equivocal", "Indeterminate", "[Not Available]", NA, ""
  ))

  expect_identical(
    as.character(observed),
    c("Negative", "Positive", "Ambiguous", "Ambiguous", "Missing", "Missing", "Missing")
  )
})

test_that("tumor classification requires three negative receptors for TNBC", {
  observed <- classify_tumor(
    er = c("Negative", "Negative", "Missing", "Positive"),
    pr = c("Negative", "Negative", "Negative", "Negative"),
    her2 = c("Negative", "Ambiguous", "Negative", "Negative")
  )

  expect_identical(
    observed,
    c("TNBC", "Unknown_tumor", "Unknown_tumor", "Other_receptor_defined")
  )
})

test_that("normal tissue stays normal regardless of tumor receptor fields", {
  metadata <- data.frame(
    file_id = c("n1", "n2", "t1"),
    file_name = c("n1.tsv", "n2.tsv", "t1.tsv"),
    path = c("n1.tsv", "n2.tsv", "t1.tsv"),
    sample_id = c("TCGA-AA-0001-11A", "TCGA-AA-0002-11A", "TCGA-AA-0003-01A"),
    patient_id = c("TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003"),
    sample_type = c("Solid Tissue Normal", "Solid Tissue Normal", "Primary Tumor")
  )
  clinical <- data.frame(
    patient_id = metadata$patient_id,
    er_status = c("Positive", NA, "Negative"),
    pr_status = c("Positive", NA, "Negative"),
    her2_status = c("Positive", NA, "Negative")
  )

  manifest <- build_sample_manifest(metadata, clinical)

  expect_identical(manifest$analysis_group, c("Normal", "Normal", "TNBC"))
  expect_true(all(manifest$eligible))
})

test_that("ambiguous tumors are excluded with explicit receptor evidence", {
  metadata <- data.frame(
    file_id = "t1",
    file_name = "t1.tsv",
    path = "t1.tsv",
    sample_id = "TCGA-AA-0001-01A",
    patient_id = "TCGA-AA-0001",
    sample_type = "Primary Tumor"
  )
  clinical <- data.frame(
    patient_id = "TCGA-AA-0001",
    er_status = "Negative",
    pr_status = "Negative",
    her2_status = "Equivocal"
  )

  manifest <- build_sample_manifest(metadata, clinical)

  expect_identical(manifest$analysis_group, "Unknown_tumor")
  expect_false(manifest$eligible)
  expect_match(manifest$exclusion_reason, "ambiguous_or_missing_receptor")
  expect_identical(as.character(manifest$her2_status), "Ambiguous")
})

test_that("deduplication is deterministic and pairing is calculated afterward", {
  manifest <- data.frame(
    patient_id = c("TCGA-AA-0001", "TCGA-AA-0001", "TCGA-AA-0001", "TCGA-AA-0002"),
    sample_id = c(
      "TCGA-AA-0001-01B", "TCGA-AA-0001-11A", "TCGA-AA-0001-01A", "TCGA-AA-0002-01A"
    ),
    analysis_group = c("TNBC", "Normal", "TNBC", "TNBC"),
    eligible = TRUE,
    exclusion_reason = NA_character_,
    stringsAsFactors = FALSE
  )

  selected_a <- select_one_per_patient_group(manifest)
  selected_b <- select_one_per_patient_group(manifest[c(4, 3, 1, 2), ])
  chosen_a <- sort(selected_a$sample_id[selected_a$selected_for_analysis])
  chosen_b <- sort(selected_b$sample_id[selected_b$selected_for_analysis])

  expect_identical(chosen_a, chosen_b)
  expect_true("TCGA-AA-0001-01A" %in% chosen_a)
  expect_false("TCGA-AA-0001-01B" %in% chosen_a)
  expect_true(all(selected_a$paired_tnbc_normal[selected_a$patient_id == "TCGA-AA-0001" &
                                                  selected_a$selected_for_analysis]))
  expect_false(selected_a$paired_tnbc_normal[
    selected_a$patient_id == "TCGA-AA-0002" & selected_a$selected_for_analysis
  ])
})

test_that("BCR clinical parser removes descriptive rows and standardizes fields", {
  path <- tempfile(fileext = ".txt")
  writeLines(c(
    paste(c(
      "bcr_patient_barcode", "er_status_by_ihc", "pr_status_by_ihc",
      "her2_status_by_ihc", "age_at_diagnosis", "race"
    ), collapse = "\t"),
    paste(c(
      "bcr_patient_barcode", "breast_carcinoma_estrogen_receptor_status",
      "breast_carcinoma_progesterone_receptor_status",
      "lab_proc_her2_neu_immunohistochemistry_receptor_status", "age", "race"
    ), collapse = "\t"),
    paste(rep("CDE_ID:", 6), collapse = "\t"),
    paste(c("TCGA-AA-0001", "Negative", "Negative", "Negative", "50", "WHITE"), collapse = "\t")
  ), path)

  clinical <- read_bcr_patient(path)

  expect_identical(clinical$patient_id, "TCGA-AA-0001")
  expect_identical(clinical$er_status, "Negative")
  expect_identical(clinical$age_at_diagnosis, 50)
  expect_identical(clinical$race, "WHITE")
})
