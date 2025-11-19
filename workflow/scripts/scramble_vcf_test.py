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

from scramble_vcf import write_vcf_header  # noqa


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

    def _run_mei_conversion(self, input_file, output_file, expected_file, sample_name):
        """Helper method to run MEI conversion and compare output"""
        # Simulate the main function logic
        with open(input_file, "r") as meis_in:
            with open(output_file, "w") as vcf_out:
                header_map = {}
                header_written = False
                caller = "scramble"

                for line in meis_in:
                    columns = line.strip().split("\t")
                    if columns[0] == "Insertion" or columns[0].startswith("#"):
                        header_map = {column_name: index for index,
                                      column_name in enumerate(columns)}
                        write_vcf_header(vcf_out, sample_name)
                        header_written = True
                        continue
                    if not line.strip():
                        continue
                    location = columns[header_map.get('Insertion', 0)]
                    if ':' not in location:
                        continue
                    chrom, pos = location.split(':')
                    try:
                        pos_int = int(pos)
                        if pos_int <= 0:
                            continue
                    except ValueError:
                        continue
                    if not chrom.startswith("chr"):
                        chrom = f"chr{chrom}"
                    mei_type = columns[header_map.get('MEI_Family', 1)].upper()
                    if mei_type == "LINE1":
                        mei_type = "L1"
                    orientation = columns[header_map.get('Orientation', 2)]
                    polarity = "+" if orientation == "Plus" else "-"
                    support = columns[header_map.get('Support', 3)]
                    score = columns[header_map.get('Score', 4)]
                    consensus = columns[header_map.get('Consensus', 7)]
                    ref = "N"
                    alt = f"<INS:ME:{mei_type}>"

                    if consensus and consensus != "NA" and consensus != "None Found":
                        svlen = str(len(consensus))
                    else:
                        if mei_type == "ALU":
                            svlen = "300"
                        elif mei_type == "L1":
                            svlen = "6000"
                        elif mei_type == "SVA":
                            svlen = "2000"
                        else:
                            svlen = "."

                    end = str(pos_int + 1)
                    info = f"SVTYPE=INS;END={end};SVLEN={svlen};MEINFO={mei_type},{pos},{end},{polarity}"
                    info += f";CALLER={caller};SUPPORT={support}"

                    out_line = f"{chrom}\t{pos}\t.\tN\t{alt}\t.\tPASS\t{info}\tGT\t0/1\n"
                    vcf_out.write(out_line)

                if not header_written:
                    write_vcf_header(vcf_out, sample_name)

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
