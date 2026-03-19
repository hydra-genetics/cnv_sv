#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Julia Höglund"
__copyright__ = "Copyright 2025, Julia Höglund"
__email__ = "julia.hoglund@scilifelab.uu.se"
__license__ = "GPL-3"


import sys
import os
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.abspath(os.path.join(TEST_DIR, "../../workflow/scripts"))
sys.path.insert(0, SCRIPT_DIR)

from scramble_vcf import write_vcf_header, process_meis_to_vcf  # noqa: E402


class TestWriteVcfHeader(unittest.TestCase):
    def test_write_vcf_header(self):
        sample_name = "TEST_SAMPLE_N"
        expected_file = "workflow/scripts/.tests/scramble_vcf.writeVcfHeader.expected.vcf"
        output_file = "workflow/scripts/.tests/scramble_vcf.writeVcfHeader.actual.vcf"

        with open(output_file, "w") as vcf_out:
            write_vcf_header(vcf_out, sample_name)

        # Compare with expected output (skip date line)
        with open(expected_file, "r") as expected:
            with open(output_file, "r") as actual:
                expected_lines = expected.readlines()
                actual_lines = actual.readlines()

                # Check fileformat line
                self.assertEqual(
                    expected_lines[0],
                    actual_lines[0],
                    "fileformat line differs"
                )

                # Compare remaining lines (skip date)
                self.assertEqual(
                    expected_lines[2:],
                    actual_lines[2:],
                    "VCF header content differs"
                )

        # Cleanup
        os.remove(output_file)


class TestMeiParsing(unittest.TestCase):
    def test_mei_conversion_alu(self):
        """Test conversion of ALU insertions from MEIs.txt to VCF"""
        input_file = "workflow/scripts/.tests/scramble_vcf.meiParsing.ALU.input.txt"
        expected_file = "workflow/scripts/.tests/scramble_vcf.meiParsing.ALU.expected.vcf"
        output_file = input_file.replace(".input.txt", ".actual.vcf")
        sample_name = "TEST_ALU_N"

        self._run_mei_conversion(input_file, output_file, expected_file, sample_name)

    def test_mei_conversion_empty(self):
        """Test conversion of empty MEIs.txt to VCF with header only"""
        input_file = "workflow/scripts/.tests/scramble_vcf.meiParsing.empty.input.txt"
        expected_file = "workflow/scripts/.tests/scramble_vcf.meiParsing.empty.expected.vcf"
        output_file = input_file.replace(".input.txt", ".actual.vcf")
        sample_name = "TEST_EMPTY_N"

        self._run_mei_conversion(input_file, output_file, expected_file, sample_name)

    def test_mei_conversion_line1(self):
        """Test conversion of LINE1 insertions from MEIs.txt to VCF"""
        input_file = "workflow/scripts/.tests/scramble_vcf.meiParsing.LINE1.input.txt"
        expected_file = "workflow/scripts/.tests/scramble_vcf.meiParsing.LINE1.expected.vcf"
        output_file = input_file.replace(".input.txt", ".actual.vcf")
        sample_name = "TEST_LINE1_N"

        self._run_mei_conversion(input_file, output_file, expected_file, sample_name)

    def test_mei_conversion_sva(self):
        """Test conversion of SVA insertions from MEIs.txt to VCF"""
        input_file = "workflow/scripts/.tests/scramble_vcf.meiParsing.SVA.input.txt"
        expected_file = "workflow/scripts/.tests/scramble_vcf.meiParsing.SVA.expected.vcf"
        output_file = input_file.replace(".input.txt", ".actual.vcf")
        sample_name = "TEST_SVA_N"

        self._run_mei_conversion(input_file, output_file, expected_file, sample_name)

    def test_mei_conversion_clustered(self):
        """Test that two nearby records with opposite Clipped_Side are clustered into one"""
        input_file = "workflow/scripts/.tests/scramble_vcf.meiParsing.clustered.input.txt"
        expected_file = "workflow/scripts/.tests/scramble_vcf.meiParsing.clustered.expected.vcf"
        output_file = input_file.replace(".input.txt", ".actual.vcf")
        sample_name = "TEST_CLUSTERED_N"

        self._run_mei_conversion(input_file, output_file, expected_file, sample_name)

    def _run_mei_conversion(self, input_file, output_file, expected_file, sample_name):
        """Helper method to run MEI conversion and compare output"""
        with open(output_file, "w") as vcf_out:
            process_meis_to_vcf(input_file, vcf_out, sample_name)

        # Compare output with expected (skip date line)
        with open(expected_file, "r") as expected:
            with open(output_file, "r") as actual:
                expected_lines = [line for line in expected.readlines() if not line.startswith("##fileDate")]
                actual_lines = [line for line in actual.readlines() if not line.startswith("##fileDate")]

                self.assertEqual(
                    expected_lines,
                    actual_lines,
                    "VCF output differs"
                )

        # Cleanup
        os.remove(output_file)


if __name__ == "__main__":
    unittest.main()
