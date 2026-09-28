# Provenance Record Sanitization (2026-09-25)

## What Changed

The checksum manifests written by the validated 2026-08-20 run recorded
absolute, machine-specific paths. Before public release, the absolute project
prefix was removed from these seven files so every entry is relative to the
repository root:

- `cohort_outputs.sha256`
- `de_outputs.sha256`
- `depmap_outputs.sha256`
- `enrichment_outputs.sha256`
- `input_checksums.sha256`
- `output_manifest.sha256`
- `source_checksums.sha256`

The transformation was a single fixed-prefix deletion applied with `sed`. No
recorded file hash changed except the one entry described below, and no entry
was added or removed.

## One Recomputed Entry

`output_manifest.sha256` records the hash of `source_checksums.sha256`. Removing
the path prefix changed that file's bytes, so its recorded hash was recomputed.

| File | Hash as recorded by the 2026-08-20 run | Hash after sanitization |
| --- | --- | --- |
| `source_checksums.sha256` | `6ddd7e6e7e7055f49af74c8670b0b594bb42efba53306e3d78e16c96e55177b9` | `dbc1c26680cf8f13b0b60484f0642e06c991712dc1831267277a91dd88fb3d73` |

## What Can Be Verified Publicly

From the repository root:

```bash
sha256sum -c analysis_v2/provenance/source_checksums.sha256
```

All 30 entries match, which shows the committed code, configuration, tests,
report template, design specification, and implementation plan are
byte-identical to the versions recorded by the validated run. The other
manifests describe frozen inputs and generated outputs that are not
redistributed (see `docs/DATA_SOURCES.md`); they can be verified only by a user
who obtains the same inputs and reruns the workflow.

## Known Limitation

`write_checksum_manifest()` and `write_output_manifest()` still write absolute
paths. The next full run will reintroduce machine-specific paths unless the
writers are changed to record project-relative paths. That change was
deliberately deferred so the published code remains byte-identical to the
validated run. It is tracked as finding O2 in `docs/AUDIT_REMEDIATION.md`.
