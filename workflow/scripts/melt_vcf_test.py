#!/usr/bin/env python3
import os
import sys
import pytest

from melt_vcf import fix_melt_vcf

# Add script directory to path
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TEST_DIR)


def test_fix_melt_vcf(tmp_path):
    input_vcf = os.path.join(TEST_DIR, ".tests/melt_vcf.input.vcf")
    expected_vcf = os.path.join(TEST_DIR, ".tests/melt_vcf.expected.vcf")
    output_vcf = tmp_path / "actual_output.vcf"

    # Run the function
    fix_melt_vcf(str(input_vcf), str(output_vcf))

    # Compare results
    with open(expected_vcf, "r") as f_exp, open(output_vcf, "r") as f_act:
        exp_lines = f_exp.readlines()
        act_lines = f_act.readlines()

        # Check total line count
        assert len(exp_lines) == len(act_lines), f"Line count mismatch: {len(exp_lines)} vs {len(act_lines)}"

        # Compare lines
        for i, (exp, act) in enumerate(zip(exp_lines, act_lines)):
            assert exp == act, f"Line {i+1} difference:\nExp: {exp.strip()}\nAct: {act.strip()}"
