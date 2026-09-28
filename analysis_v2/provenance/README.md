# Corrected Analysis Provenance

`gdc_expression_metadata.tsv` is created once by the explicitly authorized
metadata-only GDC lookup. It maps the UUID directories of the existing local
STAR-count files to TCGA samples. Subsequent runs load this table offline and
fail if it is absent or inconsistent with the local files.

`input_checksums.sha256` records the frozen expression, clinical, and DepMap
inputs. Generated run logs and package/session records remain local because
they describe a particular execution environment.

The committed manifests use repository-relative paths. They were sanitized
after the validated run; see `SANITIZATION.md` for the exact change and the one
recomputed hash.
