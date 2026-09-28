# TNBC Target Discovery: Restart Status

## Current State

The corrected frozen-data workflow under `analysis_v2/` is implemented through
all six stages. It preserves the 2025 exploratory scripts and outputs rather
than rewriting them. Source modules, configuration, tests, numbered stages,
checksums, figures, candidate evidence, and a self-contained scientific report
are now present. The final clean `run_all.R` release-gate execution completed
successfully on 2026-08-20.

## Corrected Cohort and Models

The authoritative cohort contains 115 unique TNBC tumors, 113 solid-tissue
normals, and 862 receptor-defined other primary tumors. Another 120 tumors were
excluded for ambiguous or missing receptor evidence; 14 duplicate patient-group
samples were excluded deterministically. Eleven patients contribute matched
TNBC-normal pairs.

| Contrast | Samples | Tested genes | FDR < 0.05 | Up | Down |
| --- | ---: | ---: | ---: | ---: | ---: |
| TNBC vs normal, age-adjusted | 115 vs 113 | 26,558 | 20,919 | 5,125 | 3,470 |
| Paired TNBC vs normal | 11 vs 11 | 25,951 | 9,692 | 2,382 | 4,021 |
| TNBC vs receptor-defined other, age-adjusted | 115 vs 862 | 25,812 | 19,264 | 2,638 | 2,247 |

Race and tissue-source site were rejected as group-confounded rather than
forced into the design. Ensembl IDs remain the unique primary keys.

## Enrichment, DepMap, and Evidence Classes

Local GO mapping covered 15,384 of 26,159 tested genes with nonmissing FDR
(58.81%); KEGG legacy covered 4,478 (17.12%). GO comes from
`org.Hs.eg.db`/`GO.db` 3.20.0 and KEGG legacy from MSigDB 2024.1.Hs. No online
enrichment endpoint was used.

DepMap integration includes 25 measured TNBC and 21 explicitly ER-positive or
HER2-positive breast comparator models, summarizes 17,916 genes, and maps
15,971 primary DE rows. The exact DepMap release remains
`unresolved_local_snapshot`. The current evidence table contains 4
`A_multistream`, 353 `B_expression_dependency`, 4,768 `C_expression_only`, and
21,433 `not_candidate` rows. The four Class A nominations are KIF2C, PFN1, PGK1,
and C19orf53. These are transparent screening classes, not validated targets.

## Validation and Reproduction

The corrected suite currently contains 53 test blocks and 143 assertions, with
no failures, errors, warnings, or skips in the final verification. Run:

```bash
module load R/4.4.2
TZ=UTC Rscript analysis_v2/tests/testthat.R
TZ=UTC Rscript analysis_v2/run_all.R
```

Primary deliverables are
`analysis_v2/results/candidate_evidence.csv` and
`analysis_v2/results/corrected_analysis_report.html`. Provenance is under
`analysis_v2/provenance/`; generated tables and plots are intentionally ignored
by Git.

## Preserved Legacy Snapshot

The original 2025 workflow compared 116 annotated TNBC samples with only 27
normals and used a single MDA-MB-231 screen plus heuristic ranking. Its tables
remain historical exploratory results. Checksums in
`provenance/baseline_2025-09-22.sha256` and
`provenance/restart_2026-08-19.sha256` protect the preserved state.

## Remaining Scientific Boundaries

- Bulk cross-sectional expression cannot establish causal dependency or efficacy.
- Only 11 matched pairs support the paired sensitivity analysis.
- TNBC is defined from clinical receptors, not a molecular subtype classifier.
- Residual site, batch, purity, and tissue-composition effects may remain.
- DepMap joins use gene symbols and the local release identifier is unresolved.
- Tractability is not assessed because no versioned curated source was supplied.
