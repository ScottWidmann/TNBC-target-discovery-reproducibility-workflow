# Audit and Remediation Record

This record lists the defects found when the 2025 exploratory workflow
(`tnbc_analysis_claude.r`) was audited on 2026-08-19, how each was corrected in
`analysis_v2/`, and which automated test enforces the correction. It also lists
findings that remain open. Every entry cites a location in this repository so a
reader can check it directly.

The legacy script and its outputs are preserved unchanged as historical
evidence. They are protected by `provenance/baseline_2025-09-22.sha256`.

## Remediated Findings

| ID | Legacy defect and location | Correction in `analysis_v2/` | Enforcing test(s) |
| --- | --- | --- | --- |
| L1 | Solid-tissue normal samples were relabeled as `Hormone_Positive` or `HER2_Positive` when their patient's tumor receptor fields were positive, because subtype assignment ran after the normal assignment (`tnbc_analysis_claude.r:236-247`). Only 27 normals remained. | Normal eligibility comes from TCGA sample type and is independent of receptor fields (`R/cohort.R`, `build_sample_manifest`). The corrected cohort has 113 normals. | `normal tissue stays normal regardless of tumor receptor fields` |
| L2 | TNBC status was assigned by exact string match with no explicit handling of equivocal, indeterminate, or missing receptor values (`tnbc_analysis_claude.r:229-234`). | Receptor values are normalized to negative, positive, ambiguous, or missing; TNBC requires three negatives; ambiguous tumors are excluded with the evidence recorded (`R/cohort.R`). 120 tumors were excluded. | `receptor normalization distinguishes negative, positive, ambiguous, and missing`; `tumor classification requires three negative receptors for TNBC`; `ambiguous tumors are excluded with explicit receptor evidence` |
| L3 | No explicit one-sample-per-patient rule was applied before modeling. | Deterministic selection of one sample per patient and group, with pairing computed afterward (`R/cohort.R`, `select_one_per_patient_group`). 14 duplicates were excluded. | `deduplication is deterministic and pairing is calculated afterward` |
| L4 | Genes were keyed by `make.unique(gene_name)`, which drops Ensembl IDs and appends synthetic suffixes to duplicate symbols (`tnbc_analysis_claude.r:279`). | Version-stripped Ensembl IDs are the unique primary key throughout; symbols are non-unique annotations (`R/expression_io.R`, `R/evidence.R`). | `STAR parser rejects duplicate version-stripped Ensembl IDs`; `gene joins retain Ensembl primary keys and duplicate symbols`; `DepMap join preserves unique Ensembl keys and duplicate gene symbols`; `evidence assembly preserves duplicate symbols and unique Ensembl keys` |
| L5 | Low-expression filtering was computed across all 1,104 samples rather than the samples in each contrast (`tnbc_analysis_claude.r:271-273`). | Filtering is computed per contrast from participating samples only (`R/differential_expression.R`, `filter_counts_for_contrast`). | `contrast filtering ignores samples outside the contrast`; `contrast filtering rejects absent samples and invalid thresholds` |
| L6 | The model was `~ subtype` with no assessment of covariates (`tnbc_analysis_claude.r:299`). | Covariates are admitted only after completeness, variation, and non-confounding checks. Age was admitted; race and tissue-source site were rejected as group-confounded (`R/differential_expression.R`, `assess_covariates`). Paired and specificity contrasts are reported separately. | `covariate gate accepts variation and rejects missingness or confounding`; `DE designs place patient before group only for paired models` |
| L7 | GO and KEGG enrichment used default backgrounds rather than the tested genes, and KEGG was queried online at run time (`tnbc_analysis_claude.r:533-543`). | The background is the set of tested genes that map to each frozen local collection; mapping losses are reported; no online endpoint is used (`R/enrichment.R`). | `ORA background contains only tested and mapped genes`; `mapping audit reports tested, mapped, and unmapped identifiers` |
| L8 | Druggability was inferred from gene-name regular expressions. For example, any symbol ending in `K` scored as a kinase and any ending in `R` as a receptor (`tnbc_analysis_claude.r:588-592`). These scores were combined with fold change and significance in a fixed-weight priority score (`tnbc_analysis_claude.r:612-616`). | Regex druggability and arithmetic scoring were removed. Evidence streams stay in separate columns and are combined only into transparent, rule-based classes (`R/evidence.R`, `assign_evidence_class`). Tractability is reported as not assessed until a versioned source is supplied. | `candidate evidence is classified without an arithmetic score`; `common essentials cannot enter the multistream class` |
| L9 | DepMap evidence came from a single MDA-MB-231 screen (`DepMapCrispr/depmap_integration_outputs/`). | A panel of 25 TNBC models is compared with 21 explicitly ER-positive or HER2-positive breast models (`R/depmap.R`). | `only explicit breast subtype labels enter comparator groups`; `dependency summaries retain panel size and selectivity direction` |
| L10 | No manifest linked the downloaded expression files to TCGA samples. | One authorized metadata-only GDC lookup was frozen with its retrieval time and is reused offline (`provenance/gdc_expression_metadata.tsv`). | `offline metadata loading never calls the query function`; `offline metadata loading fails when the frozen table is absent`; `metadata initialization writes retrieval time and refuses overwrite` |

## Open Findings

These are retained deliberately rather than hidden. Scientific limitations
reported to users of the results are also recorded in
`analysis_v2/provenance/warnings.csv` (W001 to W007).

| ID | Finding | Evidence | Status |
| --- | --- | --- | --- |
| O1 | Frozen-input checksums are recorded but not compared. The design specification says stages fail on checksum mismatch, but stage 00 rewrites `input_checksums.sha256` on every run and no code reads a manifest back. A changed input file with an unchanged UUID set would not be detected. | `analysis_v2/scripts/00_verify_inputs.R:99-102`; `docs/superpowers/specs/2026-08-19-corrected-tnbc-analysis-design.md`, "Execution and Failure Behavior" | Open. Proposed fix: compare against the committed manifest before writing and fail on any difference, with a fixture test. |
| O2 | Checksum manifests are written with absolute paths, which exposes machine-specific paths and prevents verification elsewhere. | `analysis_v2/R/io_contracts.R:75-94`; `analysis_v2/R/reporting.R:1-30` | Partly remediated. Committed records were sanitized (`analysis_v2/provenance/SANITIZATION.md`); writer change deferred to the next validated run. |
| O3 | The exact DepMap release of the local snapshot cannot be established from the files. | W001 | Open and disclosed as `unresolved_local_snapshot`. |
| O4 | Tractability is not assessed. | W002; `analysis_v2/config/tractability_evidence.csv` (empty schema) | Open by design until a versioned curated source is supplied. |
| O5 | DepMap evidence is joined to expression results by gene symbol, not Ensembl ID. | `analysis_v2/R/depmap.R`, `join_depmap_to_de` | Open and disclosed. Duplicate symbols are preserved rather than collapsed. |

## Disclosure Note

`provenance/restart_2026-08-19.sha256` is a point-in-time record. Two of its
eleven entries, `PROJECT_RESTART_STATUS.md` and `.gitignore`, were updated
after it was written to record the completed run and exclude generated
outputs. The other nine entries verify.
