# Provenance Records

`baseline_2025-09-22.sha256` records the code, DepMap inputs, TCGA result
tables, and plots that existed before the 2026-08-19 restart. The date in the
filename denotes the latest analysis activity in that baseline, not the date
the checksum file was created.

The 4.9 GB `GDCdata/` tree was not hashed file-by-file. Its existing download
manifest is represented by the `MANIFEST.txt` checksum. All large DepMap input
files used by the repaired integration are hashed directly.

The baseline hash for `DepMapCrispr/DepMapCrispr.R` is expected to differ after
the restart because that script was repaired. Other baseline entries should
remain unchanged.

Verify preserved baseline artifacts with:

```bash
awk '$2 != "DepMapCrispr/DepMapCrispr.R"' provenance/baseline_2025-09-22.sha256 | sha256sum -c -
```
