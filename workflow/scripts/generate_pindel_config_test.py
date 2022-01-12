#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Martin Rippin"
__copyright__ = "Copyright 2021, Martin Rippin"
__email__ = "martin.rippin@igp.uu.se"
__license__ = "GPL-3"


import io
import os
import unittest
from dataclasses import dataclass
from generate_pindel_config import insertSize, writeConfigFile


class TestInsertSize(unittest.TestCase):
    def test_insertSize(self):
        @dataclass
        class TestCase:
            name: str
            input: str
            expected: str

        testcases = [
                TestCase(
                    name="read metrics file successfully",
                    input=(
                        "workflow/scripts/.tests/"
                        "generate_pindel_config.insertSize.tsv"
                    ),
                    expected="200",
                ),
        ]

        for case in testcases:
            actual = insertSize(case.input)
            self.assertEqual(
                case.expected,
                actual,
                "failed test '{}': expected {}, got {}".format(
                    case.name, case.expected, actual
                ),
            )


class TestWriteConfigFile(unittest.TestCase):
    def test_writeConfigFile(self):
        @dataclass
        class TestCase:
            name: str
            output: str
            input: str
            insert_size: str
            sample_id: str
            expected: str

        testcases = [
                TestCase(
                    name="write config file successfully",
                    output=(
                        "workflow/scripts/.tests/"
                        "generate_pindel_config.writeConfigFile.actual.tsv"
                    ),
                    input="test.bam",
                    insert_size="200",
                    sample_id="test",
                    expected=(
                        "workflow/scripts/.tests/"
                        "generate_pindel_config.writeConfigFile.expected.tsv"
                    ),
                ),
        ]

        for case in testcases:
            writeConfigFile(
                case.output,
                case.input,
                case.insert_size,
                case.sample_id
            )
            self.assertListEqual(
                list(io.open(case.expected)),
                list(io.open(case.output)),
                "failed test '{}': files are different".format(
                    case.name
                ),
            )
            os.remove(case.output)
