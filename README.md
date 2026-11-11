# NSCLC metadata quality control

A Python workflow for extracting, standardizing and validating sample metadata from [GEO study GSE119004](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119004), a study of non-small cell lung cancer organoids.

The study includes patient tumour tissue, organoids and patient-derived xenografts. These sample types need consistent labels and documented identifiers before their metadata can be used reliably in downstream analysis.

## Workflow

1. Parse the compressed GEO SOFT file, retaining repeated sample attributes and database links.
2. Build a sample table and separate tables for characteristics and relations.
3. Standardize sample classes, tissue, organism and assay labels while retaining the submitted values.
4. Check required fields, accession formats, duplicates, series membership and agreement between source fields.
5. Export the standardized metadata, QC issues, missingness summary and transformation log.

Each run also records software versions and checksums so that its inputs and outputs can be traced.

## Results

The saved source contains 23 GEO sample records.

| Sample class | Records |
|---|---:|
| Organoid | 9 |
| Patient tumour | 8 |
| Patient-derived xenograft | 6 |
| **Total** | **23** |

The records passed the implemented QC rules with **0 errors and 0 warnings**. Pandera validation passed for the five schema fields. Four deliberately introduced defects were detected in a separate test copy.

The automated suite reported **28 tests passed and 6 subtests passed** under pytest 9.2.0. See [validation results](docs/validation_results.md) for the recorded environment and source checksum.

## Run locally

The workflow was tested with Python 3.14.7 and pandas 3.1.1.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt -r requirements-optional.txt

mkdir -p data/raw
curl --fail --location --output data/raw/GSE119004_family.soft.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119004/soft/GSE119004_family.soft.gz

python -m pytest -q
python metadata_qc.py
```

Outputs are written to a new timestamped directory under `outputs/`. A run with error-level QC findings returns a nonzero exit code and saves its reports for review.

Pandera and pytest are optional dependencies. The core workflow uses pandas, and the test suite also runs with `python -m unittest discover -s tests -v`. If Pandera is unavailable, the validation report records `not_run`.

## Output files

| File | Contents |
|---|---|
| `metadata_standardized.csv` | Submitted values alongside standardized metadata |
| `characteristics_long.csv` | Repeated sample characteristics |
| `relations_long.csv` | Documented database relations |
| `qc_issues.csv` | Record-level errors and warnings |
| `missingness.csv` | Missing-value counts by field |
| `transformation_log.csv` | Changes from submitted to standardized values |
| `pandera_status.json` | Schema validation result |

The run directory also contains the extracted source records, synthetic validation examples and a manifest with checksums and package versions.

## Validation decisions

The checks are specific to this study. GEO sample accessions must be unique, while repeated BioSample or SRA relations trigger review warnings. Host mouse strain is checked in the context of PDX records. Sample class comes from documented source fields.

The [data dictionary](docs/data_dictionary.md) describes the columns, missing-value conventions and individual QC rules.

## Scope

One row represents a GEO sample record. Patient-to-model relationships remain unverified. Passing metadata checks provides no genetic authentication or evidence of complete clinical annotation. The workflow covers metadata preparation and validation only.

Raw downloads and generated results are excluded from version control. The source can be downloaded from NCBI using the command above.
