# Data dictionary and validation decisions

## Unit of analysis

One row in `metadata_standardized.csv` is one GEO sample record, not necessarily one patient or sequencing run. The extractor reconciles actual accessions and their multiplicity against the source's declared series membership.

## Source columns

| Column | Source and meaning | Handling |
|---|---|---|
| `row_number` | 1-based extraction order | Internal locator, not biological ID |
| `geo_sample_id` | `^SAMPLE` record header | Required, unique; GSM plus digits |
| `sample_accession_raw` | `!Sample_geo_accession` | Must agree with record header |
| `sample_title_raw` | Submitted sample title | Preserve; do not infer patient links |
| `source_name_raw` | Source-name field | Independent classification evidence |
| `sample_type_raw` | Sample-type characteristic | Complete source description retained |
| `tissue_raw` | Tissue characteristic | Retained |
| `organism_raw` | Organism field | Retained |
| `library_strategy_raw` | Library strategy | Retained |
| `host_mouse_strain_raw` | Host-strain characteristic | Conditional relevance for PDX |
| `biosample_url_raw` | BioSample relation | Original documented URL |
| `sra_url_raw` | SRA relation | Original documented experiment URL |

If a scalar field has multiple occurrences, the table leaves it unresolved and flags the ambiguity. All occurrences remain in `source_records.json`. Repeated characteristics and relations also remain in their long tables.

## Derived columns

| Column | Derivation and meaning |
|---|---|
| `sample_class` | Explicit source-characteristic mapping to patient_tumour, organoid or PDX |
| `source_class` | Independent source-name mapping; compare with sample_class |
| `tissue` | Whitespace/case standardization; study scope expects lung |
| `organism` | Explicit standardization to Homo sapiens when the source supports it |
| `library_strategy` | Explicit standardization to RNA-Seq when supported |
| `host_mouse_strain` | Whitespace and documented missing-value handling |
| `biosample_id` | Accession extracted from documented NCBI BioSample URL |
| `sra_experiment_id` | SRX/ERX/DRX experiment accession, not SRR/ERR/DRR run |
| `geo_series_id` | Project scope GSE119004; source study is checked |
| `source_url` | GEO sample-page URL constructed from the accession; not fetched to validate it |
| `patient_model_link_status` | not_verified_from_selected_metadata; never inferred from title digits |

## Missing-value policy

For derived fields only, recognize None/pandas missing values plus blank strings, NA, N/A, null, none, not provided and not available (case-insensitively). These are explicit project conventions, not universal rules. Raw values stay unchanged. No rows are automatically deleted and no values are imputed.

Missingness profiling and QC failure are different. Host-strain absence outside PDX records is not automatically an error.

## Validation rules

Errors violate a defined project requirement. Warnings need review and may be legitimate.

| Rule | Severity | Meaning |
|---|---|---|
| MULTIPLE_SCALAR_VALUES | Error | Multiple source occurrences in a scalar field |
| UNPARSEABLE_CHARACTERISTIC / UNPARSEABLE_RELATION | Error | Cannot safely extract a labelled value |
| MISSING_COLUMN / MISSING_REQUIRED | Error | Missing table structure or required content |
| WRONG_STUDY / NO_SERIES_SAMPLE_LIST / SERIES_SAMPLE_MISMATCH | Error | Incorrect scope or failed membership reconciliation |
| INVALID_GEO_ID / ACCESSION_DISAGREEMENT | Error | Format or within-record accession disagreement |
| DUPLICATE_IDENTIFIER | Error for GEO; warning for BioSample/SRA | Duplicate sample rows violate table design; shared relations may be legitimate |
| CATEGORY_CONFLICT | Error | Independent source fields imply different categories |
| OUTSIDE_STUDY_SCOPE | Error | Unexpected tissue/species/assay for this study-specific workflow |
| PDX_HOST_MISSING | Warning | PDX record lacks expected host context |
| HOST_CONTEXT_REVIEW | Warning | Host information in another sample class needs interpretation |
| RELATION_UNRESOLVED | Warning | Missing/unparseable documented relation |
| INVALID_ACCESSION | Error | Wrong accession format or entity level |

Passing these rules does not prove complete clinical annotation, correct patient-model matching or biological identity. There is no live accession lookup or genetic authentication.

Resources: [NCBI SOFT](https://www.ncbi.nlm.nih.gov/geo/info/soft.html), [GSE119004](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119004), [pandas missing-data guide](https://pandas.pydata.org/docs/user_guide/missing_data.html).
