"""GSE119004 metadata QC: preserve the source, standardize cautiously, report exceptions.

This is a study-specific educational workflow, not a clinical validation system.
The core pipeline uses pandas and the Python standard library only.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import re
from typing import Any
from urllib.parse import parse_qs, urlsplit

import pandas as pd

STUDY_ID = "GSE119004"
SOURCE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119004/soft/GSE119004_family.soft.gz"
PIPELINE_VERSION = "1.0.0"
CLASS_MAP = {
    "patient tumor": "patient_tumour",
    "non-small cell lung cancer (nsclc) patient tumor tissue": "patient_tumour",
    "organoid": "organoid",
    "pdx": "PDX",
    "patient-derived xenograft (pdx)": "PDX",
}
MISSING_TOKENS = {"", "na", "n/a", "null", "none", "not provided", "not available"}
ISSUE_COLUMNS = ["row_number", "geo_sample_id", "severity", "rule", "field", "observed", "message"]
CHANGE_COLUMNS = ["row_number", "geo_sample_id", "raw_field", "derived_field", "before", "after", "method"]
RAW_FIELD_MAP = {
    "sample_accession_raw": "!Sample_geo_accession",
    "sample_title_raw": "!Sample_title",
    "source_name_raw": "!Sample_source_name_ch1",
    "organism_raw": "!Sample_organism_ch1",
    "library_strategy_raw": "!Sample_library_strategy",
}


def checksum(path: Path) -> str:
    """Fingerprint contents; does not prove biological correctness."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_soft_text(text: str) -> dict[str, Any]:
    """Read record boundaries and preserve every repeated sample attribute."""
    samples: list[dict[str, Any]] = []
    series: dict[str, list[str]] = {}
    series_id = None
    current = None
    entity = None
    for number, line in enumerate(text.splitlines(), start=1):
        if line.startswith("^"):
            if current is not None:
                samples.append(current)
            current = None
            key, sep, value = line.partition(" = ")
            if not sep:
                raise ValueError(f"Malformed SOFT record header at line {number}")
            entity = key
            if key == "^SERIES":
                if series_id is not None:
                    raise ValueError("This study-specific parser expects one SERIES record.")
                series_id = value.strip()
            elif key == "^SAMPLE":
                current = {"geo_sample_id": value.strip()}
        elif current is not None and line.startswith("!Sample_"):
            key, sep, value = line.partition(" = ")
            if not sep:
                raise ValueError(f"Malformed sample metadata at line {number}")
            current.setdefault(key, []).append(value)
        elif entity == "^SERIES" and line.startswith("!Series_"):
            key, sep, value = line.partition(" = ")
            if not sep:
                raise ValueError(f"Malformed series metadata at line {number}")
            series.setdefault(key, []).append(value)
    if current is not None:
        samples.append(current)
    if not samples:
        raise ValueError("No sample records were found.")
    return {"series_id": series_id, "series": series, "samples": samples}


def read_soft(path: Path) -> dict[str, Any]:
    """Read the saved gzip file without changing it."""
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return parse_soft_text(handle.read())


def clean_text(value: Any) -> str | None:
    """Canonical whitespace and missing tokens for DERIVED values only."""
    if value is None or pd.isna(value):
        return None
    value = " ".join(str(value).split())
    return None if value.casefold() in MISSING_TOKENS else value


def classify(value: Any) -> str | None:
    text = clean_text(value)
    return CLASS_MAP.get(text.casefold()) if text is not None else None


def issue(row: int | None, sample: str | None, severity: str, rule: str,
          field: str, observed: Any, message: str) -> dict[str, Any]:
    return dict(row_number=row, geo_sample_id=sample, severity=severity, rule=rule,
                field=field, observed=observed, message=message)


def _one(values: list[str], row: int, sample: str, field: str,
         issues: list[dict[str, Any]]) -> str | None:
    """A scalar needs exactly one occurrence; preserve ambiguity in raw JSON."""
    if not values:
        return None
    if len(values) != 1:
        issues.append(issue(row, sample, "error", "MULTIPLE_SCALAR_VALUES", field,
                            json.dumps(values), "Expected one value; all occurrences retained in extracted JSON. Resolve before selecting a scalar."))
        return None
    return values[0]


def extract_tables(parsed: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Flatten selected fields; retain repeated characteristics and relations."""
    rows, characteristics, relations, problems = [], [], [], []
    for row_number, sample in enumerate(parsed["samples"], start=1):
        sid = sample["geo_sample_id"]
        row = {"row_number": row_number, "geo_sample_id": sid}
        for column, key in RAW_FIELD_MAP.items():
            row[column] = _one(sample.get(key, []), row_number, sid, key, problems)
        char_values: dict[str, list[str]] = {}
        for order, raw in enumerate(sample.get("!Sample_characteristics_ch1", []), start=1):
            name, sep, value = raw.partition(":")
            key = (clean_text(name) or "").casefold() if sep else ""
            value = value.lstrip() if sep else None
            characteristics.append(dict(row_number=row_number, geo_sample_id=sid,
                                        occurrence=order, raw_entry=raw, characteristic=key,
                                        value_raw=value))
            if not key:
                problems.append(issue(row_number, sid, "error", "UNPARSEABLE_CHARACTERISTIC",
                                      "!Sample_characteristics_ch1", raw, "Expected a named 'key: value' entry; raw entry retained."))
            else:
                char_values.setdefault(key, []).append(value)
        for column, key in [("sample_type_raw", "sample type"), ("tissue_raw", "tissue"),
                            ("host_mouse_strain_raw", "host mouse strain")]:
            row[column] = _one(char_values.get(key, []), row_number, sid, key, problems)
        rel_values: dict[str, list[str]] = {}
        for order, raw in enumerate(sample.get("!Sample_relation", []), start=1):
            name, sep, value = raw.partition(":")
            label = (clean_text(name) or "").casefold() if sep else ""
            value = value.lstrip() if sep else None
            relations.append(dict(row_number=row_number, geo_sample_id=sid, occurrence=order,
                                  raw_entry=raw, relation_type=label, url_raw=value))
            if not label:
                problems.append(issue(row_number, sid, "error", "UNPARSEABLE_RELATION",
                                      "!Sample_relation", raw, "Expected 'database: URL'; raw entry retained."))
            else:
                rel_values.setdefault(label, []).append(value)
        row["biosample_url_raw"] = _one(rel_values.get("biosample", []), row_number, sid, "BioSample relation", problems)
        row["sra_url_raw"] = _one(rel_values.get("sra", []), row_number, sid, "SRA relation", problems)
        rows.append(row)
    return (pd.DataFrame(rows), pd.DataFrame(characteristics), pd.DataFrame(relations),
            pd.DataFrame(problems, columns=ISSUE_COLUMNS))


def accession_from_url(value: Any, database: str) -> str | None:
    """Extract documented NCBI accessions only; never infer patient links."""
    text = clean_text(value)
    if text is None:
        return None
    parts = urlsplit(text)
    if parts.scheme not in {"https", "http"} or parts.hostname not in {"www.ncbi.nlm.nih.gov", "ncbi.nlm.nih.gov"}:
        return None
    if database == "biosample":
        match = re.fullmatch(r"/biosample/(SAM[NED]\d+)/?", parts.path)
        return match.group(1) if match else None
    if database == "sra" and parts.path.rstrip("/") == "/sra":
        terms = parse_qs(parts.query).get("term", [])
        if len(terms) == 1 and re.fullmatch(r"[SED]RX\d+", terms[0]):
            return terms[0]
    return None


def standardize(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Add derived columns without altering raw columns or dropping records."""
    table = raw.copy(deep=True)
    mappings = {
        "sample_class": ("sample_type_raw", classify, "Explicit study-specific sample-type lookup"),
        "source_class": ("source_name_raw", classify, "Independent source-name lookup for agreement checking"),
        "tissue": ("tissue_raw", lambda v: (clean_text(v) or "").casefold() or None, "Whitespace and case standardization"),
        "organism": ("organism_raw", lambda v: "Homo sapiens" if (clean_text(v) or "").casefold() == "homo sapiens" else clean_text(v), "Explicit species-label mapping"),
        "library_strategy": ("library_strategy_raw", lambda v: "RNA-Seq" if (clean_text(v) or "").casefold() == "rna-seq" else clean_text(v), "Explicit assay-label mapping"),
        "host_mouse_strain": ("host_mouse_strain_raw", clean_text, "Whitespace and missing-token standardization"),
        "biosample_id": ("biosample_url_raw", lambda v: accession_from_url(v, "biosample"), "Extract accession from documented NCBI BioSample URL"),
        "sra_experiment_id": ("sra_url_raw", lambda v: accession_from_url(v, "sra"), "Extract SRX/ERX/DRX experiment accession; this is NOT a run accession"),
    }
    changes = []
    for derived, (original, transform, reason) in mappings.items():
        values = []
        for _, row in raw.iterrows():
            before = row[original]
            after = transform(before)
            values.append(after)
            before_text = None if pd.isna(before) else str(before)
            if before_text != after:
                changes.append(dict(row_number=int(row["row_number"]), geo_sample_id=row["geo_sample_id"],
                                    raw_field=original, derived_field=derived, before=before_text,
                                    after=after, method=reason))
        table[derived] = pd.Series(values, index=table.index, dtype="string")
    table["geo_series_id"] = STUDY_ID
    table["source_url"] = table["geo_sample_id"].map(lambda s: f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={s}")
    table["patient_model_link_status"] = "not_verified_from_selected_metadata"
    return table, pd.DataFrame(changes, columns=CHANGE_COLUMNS)


def validate(table: pd.DataFrame, parsed: dict[str, Any], extraction_issues: pd.DataFrame | None = None) -> pd.DataFrame:
    """Return row-level errors/warnings. Zero issues means these rules passed only."""
    findings = [] if extraction_issues is None else extraction_issues.to_dict("records")
    def add(row, severity, rule, field, observed, message):
        findings.append(issue(None if row is None else int(row["row_number"]),
                              None if row is None else row["geo_sample_id"],
                              severity, rule, field, observed, message))
    required = ["geo_sample_id", "sample_accession_raw", "sample_title_raw", "source_name_raw",
                "sample_type_raw", "tissue", "organism", "library_strategy", "sample_class", "source_class"]
    needed = required + ["row_number", "host_mouse_strain", "biosample_id", "sra_experiment_id"]
    for column in needed:
        if column not in table.columns:
            add(None, "error", "MISSING_COLUMN", column, None, "Required schema column absent.")
    if any(c not in table.columns for c in needed):
        return pd.DataFrame(findings, columns=ISSUE_COLUMNS)
    if parsed["series_id"] != STUDY_ID:
        add(None, "error", "WRONG_STUDY", "series_id", parsed["series_id"], "This workflow is configured for GSE119004.")
    expected = parsed["series"].get("!Series_sample_id", [])
    actual = table["geo_sample_id"].tolist()
    if not expected:
        add(None, "error", "NO_SERIES_SAMPLE_LIST", "geo_sample_id", None, "Cannot reconcile against series membership.")
    elif Counter(actual) != Counter(expected):
        add(None, "error", "SERIES_SAMPLE_MISMATCH", "geo_sample_id", json.dumps(actual),
            "Extracted IDs/multiplicity differ from the series sample list.")
    if len(table) == 0:
        add(None, "error", "EMPTY_TABLE", "row_number", 0, "No records available for validation.")
    for _, row in table.iterrows():
        for field in required:
            if clean_text(row[field]) is None:
                add(row, "error", "MISSING_REQUIRED", field, None, "Missing, blank, ambiguous or unrecognized required value; investigate source.")
        if not re.fullmatch(r"GSM\d+", str(row["geo_sample_id"])):
            add(row, "error", "INVALID_GEO_ID", "geo_sample_id", row["geo_sample_id"], "Expected GSM followed by digits; no automatic ID repair.")
        if clean_text(row["sample_accession_raw"]) != clean_text(row["geo_sample_id"]):
            add(row, "error", "ACCESSION_DISAGREEMENT", "sample_accession_raw", row["sample_accession_raw"], "Record header differs from declared sample accession.")
        sample_class, source_class = clean_text(row["sample_class"]), clean_text(row["source_class"])
        if sample_class is not None and source_class is not None and sample_class != source_class:
            add(row, "error", "CATEGORY_CONFLICT", "sample_class", sample_class, "Sample type and source name imply different classes.")
        for field, allowed in [("tissue", "lung"), ("organism", "Homo sapiens"), ("library_strategy", "RNA-Seq")]:
            if clean_text(row[field]) is not None and row[field] != allowed:
                add(row, "error", "OUTSIDE_STUDY_SCOPE", field, row[field], f"This study-specific workflow expects {allowed}; do not generalize to all GEO datasets.")
        if sample_class == "PDX" and clean_text(row["host_mouse_strain"]) is None:
            add(row, "warning", "PDX_HOST_MISSING", "host_mouse_strain", None, "Host strain should be reviewed for a PDX record.")
        elif sample_class not in {None, "PDX"} and clean_text(row["host_mouse_strain"]) is not None:
            add(row, "warning", "HOST_CONTEXT_REVIEW", "host_mouse_strain", row["host_mouse_strain"], "Host information outside PDX records may need contextual review; not automatically wrong.")
        for field, pattern in [("biosample_id", r"SAM[NED]\d+"), ("sra_experiment_id", r"[SED]RX\d+")]:
            value = clean_text(row[field])
            if value is None:
                add(row, "warning", "RELATION_UNRESOLVED", field, None, "Missing or unparseable documented relation; no relationship inferred.")
            elif not re.fullmatch(pattern, value):
                add(row, "error", "INVALID_ACCESSION", field, value, "Accession does not match the expected identifier type.")
    for field, severity in [("geo_sample_id", "error"), ("biosample_id", "warning"), ("sra_experiment_id", "warning")]:
        values = table[field].map(clean_text)
        duplicate = values.notna() & values.duplicated(keep=False)
        for _, row in table.loc[duplicate].iterrows():
            add(row, severity, "DUPLICATE_IDENTIFIER", field, row[field],
                "Repeated GEO record is an error; shared BioSample/SRA links require review and may be valid. No rows deleted.")
    return pd.DataFrame(findings, columns=ISSUE_COLUMNS)


def missingness_report(table: pd.DataFrame) -> pd.DataFrame:
    """Keep generic missingness separate from context-aware QC failures."""
    rows = []
    for field in table.columns:
        missing = table[field].map(clean_text).isna()
        rows.append({"field": field, "missing_count": int(missing.sum()), "total_rows": len(table),
                     "missing_percent": round(float(missing.mean() * 100), 2) if len(table) else 0})
    return pd.DataFrame(rows)


def synthetic_demo(table: pd.DataFrame, parsed: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Inject errors into a COPY. These are NOT observed study defects."""
    if len(table) < 3:
        raise ValueError("The demo needs at least three records.")
    demo = table.copy(deep=True)
    a, b, c = demo.index[:3]
    demo.loc[a, "tissue"] = pd.NA
    demo.loc[b, "geo_sample_id"] = demo.loc[a, "geo_sample_id"]
    demo.loc[c, "source_class"] = "organoid" if demo.loc[c, "sample_class"] != "organoid" else "PDX"
    demo.loc[c, "sra_experiment_id"] = "SRR123456"
    manifest = pd.DataFrame([
        {"row_number": int(demo.loc[a,"row_number"]), "injected_change": "Remove tissue", "expected_rule": "MISSING_REQUIRED"},
        {"row_number": int(demo.loc[b,"row_number"]), "injected_change": "Duplicate a GEO sample identifier", "expected_rule": "DUPLICATE_IDENTIFIER"},
        {"row_number": int(demo.loc[c,"row_number"]), "injected_change": "Conflicting source/sample category", "expected_rule": "CATEGORY_CONFLICT"},
        {"row_number": int(demo.loc[c,"row_number"]), "injected_change": "Use an SRR run ID where SRX experiment is expected", "expected_rule": "INVALID_ACCESSION"},
    ])
    return demo, validate(demo, parsed), manifest


def pandera_check(table: pd.DataFrame) -> dict[str, Any]:
    """Optional schema validation; explicitly skipped when unavailable."""
    try:
        import pandera.pandas as pa
    except ImportError as exc:
        return {"status": "not_run", "reason": "Optional Pandera dependency unavailable; no packages were installed.", "detail": str(exc)}
    schema = pa.DataFrameSchema({
        "geo_sample_id": pa.Column(str, pa.Check.str_matches(r"^GSM\d+$"), unique=True, nullable=False),
        "sample_class": pa.Column(str, pa.Check.isin(["patient_tumour", "organoid", "PDX"]), nullable=False),
        "tissue": pa.Column(str, pa.Check.isin(["lung"]), nullable=False),
        "organism": pa.Column(str, pa.Check.isin(["Homo sapiens"]), nullable=False),
        "library_strategy": pa.Column(str, pa.Check.isin(["RNA-Seq"]), nullable=False),
    }, strict=False, coerce=False)
    try:
        schema.validate(table, lazy=True)
        return {"status": "passed", "scope": "Selected schema fields; does not replace custom semantic QC."}
    except pa.errors.SchemaErrors as exc:
        return {"status": "failed", "failures": json.loads(exc.failure_cases.to_json(orient="records"))}


def run_pipeline(project_dir: Path) -> dict[str, Any]:
    """Run from the preserved source and write a NEW timestamped output directory."""
    root = project_dir.resolve()
    source = root / "data" / "raw" / "GSE119004_family.soft.gz"
    before_hash = checksum(source)
    parsed = read_soft(source)
    raw, characteristics, relations, extraction_issues = extract_tables(parsed)
    table, changes = standardize(raw)
    qc = validate(table, parsed, extraction_issues)
    missing = missingness_report(table)
    demo, demo_qc, injection = synthetic_demo(table, parsed)
    pandera_result = pandera_check(table)
    now = datetime.now(timezone.utc)
    out = root / "outputs" / now.strftime("%Y%m%dT%H%M%S_%fZ")
    out.mkdir(parents=True, exist_ok=False)
    frames = {
        "metadata_extracted.csv": raw,
        "metadata_standardized.csv": table,
        "characteristics_long.csv": characteristics,
        "relations_long.csv": relations,
        "transformation_log.csv": changes,
        "qc_issues.csv": qc,
        "missingness.csv": missing,
        "synthetic_demo_metadata.csv": demo,
        "synthetic_demo_qc_issues.csv": demo_qc,
        "synthetic_injections.csv": injection,
    }
    for name, frame in frames.items():
        frame.to_csv(out / name, index=False, na_rep="")
    (out / "source_records.json").write_text(json.dumps(parsed, indent=2, ensure_ascii=False), encoding="utf-8")
    (out / "pandera_status.json").write_text(json.dumps(pandera_result, indent=2), encoding="utf-8")
    errors = int((qc["severity"] == "error").sum())
    warnings = int((qc["severity"] == "warning").sum())
    summary = {
        "sample_records": len(table), "unique_geo_ids": int(table["geo_sample_id"].nunique()),
        "sample_class_counts": {str(k): int(v) for k,v in table["sample_class"].value_counts(dropna=False).items()},
        "qc_errors": errors, "qc_warnings": warnings,
        "status": "FAIL" if errors else ("PASS_WITH_WARNINGS" if warnings else "PASS_SELECTED_RULES"),
        "synthetic_injections": len(injection),
        "synthetic_expected_rules_detected": sorted(set(injection["expected_rule"]) & set(demo_qc["rule"])),
        "pandera_status": pandera_result["status"],
    }
    versions = {}
    for name in ("pandas", "pandera", "pytest", "ipykernel", "certifi"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "not installed"
    manifest = {
        "study_id": STUDY_ID, "source_url": SOURCE_URL,
        "source_relative_path": str(source.relative_to(root)), "source_sha256": before_hash,
        "source_bytes": source.stat().st_size, "original_download_time_utc": None,
        "download_time_note": "The original notebook did not record a download timestamp. This run timestamp is not substituted for it.",
        "pipeline_run_time_utc": now.isoformat(), "pipeline_version": PIPELINE_VERSION,
        "python_version": platform.python_version(), "package_versions": versions,
        "pipeline_sha256": checksum(Path(__file__)), "summary": summary,
        "limits": ["Metadata checks do not verify biological identity.",
                   "No patient/model relationships inferred from sample titles.",
                   "SRX accessions identify experiments, not SRR sequencing runs.",
                   "Public source labels alone do not establish direct patient versus xenograft origin of organoids.",
                   "No raw-read QC, expression normalization, clinical data integration or machine learning performed.",
                   "Synthetic defects are educational tests, not findings about the real study."],
    }
    manifest["output_sha256"] = {p.name: checksum(p) for p in sorted(out.iterdir()) if p.is_file()}
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    counts = ", ".join(f"{v} {k}" for k,v in summary["sample_class_counts"].items())
    report = f"""# GSE119004 metadata QC report

Run (UTC): {now.isoformat()}

## Results from the saved source
- {len(table)} records; {summary['unique_geo_ids']} unique GEO sample identifiers.
- Sample categories: {counts}.
- Selected-rule status: **{summary['status']}**; {errors} errors; {warnings} warnings.
- {len(characteristics)} characteristic entries and {len(relations)} relation entries retained.
- {len(changes)} logged transformations/derived-field extractions. These are not errors corrected.
- Optional Pandera: {pandera_result['status']}.

## Interpretation
A pass means the implemented metadata rules passed. It does not establish biological sample identity,
validate the experiment, confirm patient-model matching or guarantee complete clinical metadata.
A blank host-mouse-strain field for a non-PDX record is not automatically a quality failure.
Review `qc_issues.csv` for exceptions; rows are never silently deleted or imputed.
`source_records.json` preserves repeated fields. `transformation_log.csv` traces derived labels.

## Test demonstration (separate from study results)
Four defects were deliberately introduced in a copy. Expected rules detected:
{', '.join(summary['synthetic_expected_rules_detected'])}.
See `synthetic_injections.csv` and `synthetic_demo_qc_issues.csv`. Additional related flags may occur.
Do not report these artificial failures as errors found in GSE119004.

## Scope and limitations
Patient-model matching remains unverified. Organoids remain labelled `organoid`, not automatically PDO.
SRA links are experiment accessions; no run-level FASTQ files were downloaded or analysed.
No differential expression, TPM normalization, read-depth assessment, tumour identity verification or AI
model was run. Original download time was not logged; the manifest records the actual pipeline run time.

## Provenance
Source: {SOURCE_URL}
SHA-256: `{before_hash}`
Python: {platform.python_version()}; pandas: {pd.__version__}
A new output directory is created on every run. Original source data remain unchanged.

## Interview takeaway
“I standardize documented metadata using explicit mappings, preserve the original labels, and
report exceptions instead of guessing. I test validation rules on deliberately faulty copies,
and keep those tests separate from real study findings.”
"""
    (out / "qc_summary.md").write_text(report, encoding="utf-8")
    if checksum(source) != before_hash:
        raise RuntimeError("Source file changed during execution.")
    return {"output_dir": out, "summary": summary, "table": table, "qc": qc,
            "manifest": manifest, "changes": changes}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    result = run_pipeline(args.project_dir)
    print(json.dumps(result["summary"], indent=2))
    print("Outputs:", result["output_dir"])
    raise SystemExit(1 if result["summary"]["qc_errors"] else 0)
