# Validation results

The workflow was run against the saved GSE119004 SOFT snapshot on 10 November 2026 (America/Toronto).

| Check | Result |
|---|---|
| Extracted sample records | 23 |
| Unique GEO sample accessions | 23 |
| Sample classes | 9 organoids, 8 patient tumours, 6 PDXs |
| Custom QC errors | 0 |
| Custom QC warnings | 0 |
| Pandera schema validation | Passed |
| pytest | 28 tests passed, 6 subtests passed |
| Synthetic defect detection | All four expected rule types detected |

The synthetic cases cover conflicting categories, duplicate identifiers, invalid accessions and missing required values. They are deliberately modified copies used to check validator behaviour.

## Environment

| Component | Version |
|---|---|
| Python | 3.14.7 |
| pandas | 3.1.1 |
| Pandera | 0.36.0 |
| pytest | 9.2.0 |

## Source

Study accession — GSE119004

Source file — `GSE119004_family.soft.gz`

SHA-256 — `66fb8deec1ccb53d4ad76c9cf1fb8716719be3ab3ec3d10f0af76462875524a1`

The original download timestamp was unavailable. The manifest records the validation run time separately.

## Interpretation

These results apply to the implemented metadata checks and this source snapshot. Biological sample identity and patient-to-model relationships require additional evidence.
