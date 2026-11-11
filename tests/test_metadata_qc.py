"""Synthetic test fixtures only. Nothing here claims errors exist in GEO."""
from copy import deepcopy
import unittest
import pandas as pd
from metadata_qc import (accession_from_url, clean_text, extract_tables, parse_soft_text,
                         standardize, synthetic_demo, validate)


def fixture_text() -> str:
    """Small invented records, distinct from the real source."""
    header = "^SERIES = GSE119004\n" + "".join(f"!Series_sample_id = GSM{i}\n" for i in (1,2,3))
    text = header
    for i, source, sample_type in [(1,"patient tumor","non-small cell lung cancer (NSCLC) patient tumor tissue"),
                                  (2,"organoid","organoid"), (3,"PDX","patient-derived xenograft (PDX)")]:
        text += (f"^SAMPLE = GSM{i}\n!Sample_geo_accession = GSM{i}\n"
                 f"!Sample_title = SYNTHETIC_{i}\n!Sample_source_name_ch1 = {source}\n"
                 "!Sample_organism_ch1 = Homo sapiens\n"
                 f"!Sample_characteristics_ch1 = sample type: {sample_type}\n"
                 "!Sample_characteristics_ch1 = tissue: lung\n"
                 "!Sample_library_strategy = RNA-Seq\n"
                 f"!Sample_relation = BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN{i}\n"
                 f"!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX{i}\n")
        if i == 3:
            text += "!Sample_characteristics_ch1 = host mouse strain: NOD/SCID\n"
    return text


class MetadataQCTests(unittest.TestCase):
    def setUp(self):
        self.parsed = parse_soft_text(fixture_text())
        self.raw, self.characteristics, self.relations, self.extraction_issues = extract_tables(self.parsed)
        self.table, self.changes = standardize(self.raw)

    def codes(self, table=None, parsed=None):
        frame = validate(self.table if table is None else table, self.parsed if parsed is None else parsed)
        return set(frame["rule"])

    def test_valid_fixture_has_no_issues(self):
        self.assertEqual(self.codes(), set())

    def test_parser_saves_last_sample_at_eof(self):
        self.assertEqual(len(self.parsed["samples"]), 3)
        self.assertEqual(self.parsed["samples"][-1]["geo_sample_id"], "GSM3")

    def test_parser_flushes_before_non_sample_entity(self):
        parsed = parse_soft_text(fixture_text() + "^PLATFORM = GPL1\n!Platform_title = not a sample\n")
        self.assertEqual(len(parsed["samples"]), 3)
        self.assertNotIn("!Platform_title", parsed["samples"][-1])

    def test_repeated_characteristics_and_relations_retained(self):
        self.assertEqual(len(self.characteristics), 7)
        self.assertEqual(len(self.relations), 6)
        self.assertEqual(len(self.parsed["samples"][0]["!Sample_relation"]), 2)

    def test_no_samples_fails_explicitly(self):
        with self.assertRaisesRegex(ValueError, "No sample"):
            parse_soft_text("^SERIES = GSE119004\n")

    def test_malformed_attribute_fails_explicitly(self):
        with self.assertRaisesRegex(ValueError, "Malformed sample"):
            parse_soft_text(fixture_text() + "!Sample_bad_delimiter\n")

    def test_characteristic_preserves_internal_colon(self):
        parsed = parse_soft_text(fixture_text() + "!Sample_characteristics_ch1 = note: a: b\n")
        _, chars, _, _ = extract_tables(parsed)
        self.assertEqual(chars.iloc[-1]["value_raw"], "a: b")

    def test_multiple_scalar_values_not_silently_chosen(self):
        parsed = deepcopy(self.parsed)
        parsed["samples"][0]["!Sample_title"].append("CONFLICT")
        raw, _, _, problems = extract_tables(parsed)
        self.assertTrue(pd.isna(raw.loc[0, "sample_title_raw"]))
        self.assertIn("MULTIPLE_SCALAR_VALUES", set(problems["rule"]))
        self.assertEqual(len(parsed["samples"][0]["!Sample_title"]), 2)

    def test_missing_characteristic_delimiter_reported(self):
        parsed = deepcopy(self.parsed)
        parsed["samples"][0]["!Sample_characteristics_ch1"].append("malformed")
        _, _, _, problems = extract_tables(parsed)
        self.assertIn("UNPARSEABLE_CHARACTERISTIC", set(problems["rule"]))

    def test_standardization_preserves_raw_columns(self):
        snapshot = self.raw.copy(deep=True)
        clean, _ = standardize(self.raw)
        pd.testing.assert_frame_equal(self.raw, snapshot)
        pd.testing.assert_frame_equal(clean[self.raw.columns], self.raw)
        self.assertEqual(len(clean), len(self.raw))

    def test_known_labels_standardize_explicitly(self):
        raw = self.raw.copy()
        raw.loc[0, "source_name_raw"] = "  Patient   TUMOR "
        clean, log = standardize(raw)
        self.assertEqual(clean.loc[0, "source_class"], "patient_tumour")
        self.assertEqual(clean.loc[0, "source_name_raw"], "  Patient   TUMOR ")
        self.assertTrue(len(log) > 0)

    def test_unrecognized_label_not_guessed(self):
        raw = self.raw.copy()
        raw.loc[0, "sample_type_raw"] = "something vaguely tumour-like"
        clean, _ = standardize(raw)
        self.assertTrue(pd.isna(clean.loc[0, "sample_class"]))
        self.assertIn("MISSING_REQUIRED", self.codes(clean))

    def test_blank_strings_and_na_tokens_recognized(self):
        for value in [None, pd.NA, "", "  ", "N/A", "not provided"]:
            with self.subTest(value=str(value)):
                self.assertIsNone(clean_text(value))
        self.assertEqual(clean_text(" NOD/SCID "), "NOD/SCID")

    def test_required_field_missing(self):
        self.table.loc[0, "tissue"] = " "
        self.assertIn("MISSING_REQUIRED", self.codes())

    def test_duplicate_geo_ids_detected_without_row_deletion(self):
        self.table.loc[1, "geo_sample_id"] = "GSM1"
        issues = validate(self.table, self.parsed)
        self.assertEqual(len(self.table), 3)
        self.assertEqual((issues["rule"] == "DUPLICATE_IDENTIFIER").sum(), 2)

    def test_series_membership_checks_ids_not_just_count(self):
        self.table.loc[0, "geo_sample_id"] = "GSM999"
        self.assertIn("SERIES_SAMPLE_MISMATCH", self.codes())

    def test_record_accession_mismatch(self):
        self.table.loc[0, "sample_accession_raw"] = "GSM999"
        self.assertIn("ACCESSION_DISAGREEMENT", self.codes())

    def test_category_conflict(self):
        self.table.loc[0, "source_class"] = "organoid"
        self.assertIn("CATEGORY_CONFLICT", self.codes())

    def test_pdx_requires_contextual_host_review(self):
        self.table.loc[2, "host_mouse_strain"] = pd.NA
        issues = validate(self.table, self.parsed)
        result = issues.loc[issues["rule"] == "PDX_HOST_MISSING"]
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["severity"], "warning")

    def test_non_pdx_missing_host_not_flagged(self):
        self.assertTrue(pd.isna(self.table.loc[0, "host_mouse_strain"]))
        self.assertNotIn("PDX_HOST_MISSING", self.codes())

    def test_valid_experiment_accession_extraction(self):
        self.assertEqual(accession_from_url("https://www.ncbi.nlm.nih.gov/sra?term=SRX123", "sra"), "SRX123")

    def test_run_accession_not_accepted_as_experiment(self):
        self.assertIsNone(accession_from_url("https://www.ncbi.nlm.nih.gov/sra?term=SRR123", "sra"))
        self.table.loc[0, "sra_experiment_id"] = "SRR123"
        self.assertIn("INVALID_ACCESSION", self.codes())

    def test_non_ncbi_relation_not_trusted(self):
        self.assertIsNone(accession_from_url("https://example.com/biosample/SAMN1", "biosample"))
        self.assertIsNone(accession_from_url("https://www.ncbi.nlm.nih.gov.evil.test/sra?term=SRX1", "sra"))

    def test_shared_biosample_is_review_warning_not_automatic_error(self):
        self.table.loc[1, "biosample_id"] = self.table.loc[0, "biosample_id"]
        issues = validate(self.table, self.parsed)
        results = issues.loc[(issues["rule"] == "DUPLICATE_IDENTIFIER") & (issues["field"] == "biosample_id")]
        self.assertEqual(set(results["severity"]), {"warning"})

    def test_study_scope_mismatch(self):
        self.table.loc[0, "organism"] = "Mus musculus"
        self.assertIn("OUTSIDE_STUDY_SCOPE", self.codes())

    def test_missing_column_reports_without_crash(self):
        table = self.table.drop(columns=["tissue"])
        self.assertIn("MISSING_COLUMN", self.codes(table))

    def test_no_patient_links_inferred(self):
        self.assertEqual(set(self.table["patient_model_link_status"]), {"not_verified_from_selected_metadata"})
        self.assertNotIn("patient_id", self.table.columns)
        self.assertEqual(self.table.loc[1, "sample_class"], "organoid")

    def test_synthetic_errors_detected_real_table_unchanged(self):
        before = self.table.copy(deep=True)
        _, findings, changes = synthetic_demo(self.table, self.parsed)
        self.assertTrue(set(changes["expected_rule"]).issubset(set(findings["rule"])))
        pd.testing.assert_frame_equal(before, self.table)
        self.assertEqual(self.codes(), set())


if __name__ == "__main__":
    unittest.main(verbosity=2)
