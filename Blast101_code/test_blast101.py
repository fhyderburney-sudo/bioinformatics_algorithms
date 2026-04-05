import os
import tempfile
import unittest

import blast101_cli
import process_fasta_file as pff
import smith_waterman_p as SW
import calc_bit_and_evalues as cbe


class TestCLIValidation(unittest.TestCase):
    """Tests for simple command-line input validation."""

    def test_valid_query_sequence(self):
        seq = "PWNAAPLHNFGEDFLQPYVQLQQNFSASDLEVNLEATRESHAHFSTPQALELFLNYSVTP"
        valid = "GAVLITSMCPFYWHKRDENQ"
        result = blast101_cli.validate_query_sequence(seq, valid)
        self.assertEqual(result, seq.upper())

    def test_empty_query_sequence(self):
        valid = "GAVLITSMCPFYWHKRDENQ"
        with self.assertRaises(ValueError):
            blast101_cli.validate_query_sequence("", valid)

    def test_invalid_query_sequence(self):
        valid = "GAVLITSMCPFYWHKRDENQ"
        with self.assertRaises(ValueError):
            blast101_cli.validate_query_sequence("ABCD123", valid)

    def test_dna_query_sequence_rejected(self):
        valid = "GAVLITSMCPFYWHKRDENQ"
        with self.assertRaises(ValueError):
            blast101_cli.validate_query_sequence("ACGTACGT", valid)

    def test_valid_database_file(self):
        with tempfile.NamedTemporaryFile("w", delete=False, suffix=".fasta") as tmp:
            tmp.write(">seq1\nACDEFG\n")
            tmp_path = tmp.name

        try:
            result = blast101_cli.validate_database_file(tmp_path)
            self.assertEqual(result, tmp_path)
        finally:
            os.remove(tmp_path)

    def test_missing_database_file(self):
        with self.assertRaises(ValueError):
            blast101_cli.validate_database_file("not_a_real_file.fasta")


class TestFastaProcessing(unittest.TestCase):
    """Tests for FASTA parsing/processing using small synthetic files."""

    @staticmethod
    def score_by_length(seq):
        return len(seq)

    def test_single_fasta_record(self):
        with tempfile.NamedTemporaryFile("w", delete=False, suffix=".fasta") as tmp:
            tmp.write(">seq1\nACDEFG\n")
            tmp_path = tmp.name

        try:
            res = pff.process_fasta_file(tmp_path, self.score_by_length, max_scores=10)
            self.assertEqual(len(res), 1)
            self.assertEqual(res[0][0], ">seq1")
            self.assertEqual(res[0][1], "ACDEFG")
            self.assertEqual(res[0][2], 6)
        finally:
            os.remove(tmp_path)

    def test_multiline_fasta_record(self):
        with tempfile.NamedTemporaryFile("w", delete=False, suffix=".fasta") as tmp:
            tmp.write(">seq1\nACD\nEFG\n")
            tmp_path = tmp.name

        try:
            res = pff.process_fasta_file(tmp_path, self.score_by_length, max_scores=10)
            self.assertEqual(len(res), 1)
            self.assertEqual(res[0][1], "ACDEFG")
            self.assertEqual(res[0][2], 6)
        finally:
            os.remove(tmp_path)

    def test_multiple_fasta_records(self):
        with tempfile.NamedTemporaryFile("w", delete=False, suffix=".fasta") as tmp:
            tmp.write(">seq1\nAAAA\n>seq2\nCCCCCC\n")
            tmp_path = tmp.name

        try:
            res = pff.process_fasta_file(tmp_path, self.score_by_length, max_scores=10)
            self.assertEqual(len(res), 2)

            headers = [r[0] for r in res]
            scores = [r[2] for r in res]

            self.assertIn(">seq1", headers)
            self.assertIn(">seq2", headers)
            self.assertIn(4, scores)
            self.assertIn(6, scores)
        finally:
            os.remove(tmp_path)

    def test_empty_fasta_file(self):
        with tempfile.NamedTemporaryFile("w", delete=False, suffix=".fasta") as tmp:
            tmp_path = tmp.name

        try:
            res = pff.process_fasta_file(tmp_path, self.score_by_length, max_scores=10)
            self.assertEqual(res, [])
        finally:
            os.remove(tmp_path)

    def test_missing_fasta_file(self):
        with self.assertRaises(FileNotFoundError):
            pff.process_fasta_file("missing_file.fasta", self.score_by_length, max_scores=10)


class TestSmithWaterman(unittest.TestCase):
    """Simple Smith-Waterman tests with known expected behaviour."""

    def test_identical_sequences_score_positive(self):
        score = SW.perform_smith_waterman("AAAA", "AAAA", False, False)
        self.assertGreater(score, 0)

    def test_identical_scores_higher_than_mismatch(self):
        score_same = SW.perform_smith_waterman("AAAA", "AAAA", False, False)
        score_diff = SW.perform_smith_waterman("AAAA", "RRRR", False, False)
        self.assertGreater(score_same, score_diff)

    def test_perfect_match_ground_truth(self):
        query = "PWNAAPLHNFGEDFLQPYVQLQQNFSASDLEVNLEATRESHAHFSTPQALELFLNYSVTP"
        score = SW.perform_smith_waterman(query, query, False, False)
        self.assertEqual(score, 320)


class TestBitAndEvalues(unittest.TestCase):
    """Simple tests for helper functions using lecture-style values."""

    def setUp(self):
        self.k = 6.23037
        self.scale = 29.01160

    def test_bit_score_positive(self):
        score = cbe.get_bit_score(300, self.k, self.scale)
        self.assertGreater(score, 0)

    def test_expect_string_not_empty(self):
        expect = cbe.get_expect_s(300, self.k, self.scale)
        self.assertTrue(len(expect) > 0)

    def test_higher_score_gives_lower_expect(self):
        low_expect = float(cbe.get_expect_s(40, self.k, self.scale))
        high_expect = float(cbe.get_expect_s(300, self.k, self.scale))
        self.assertLess(high_expect, low_expect)


if __name__ == "__main__":
    unittest.main()