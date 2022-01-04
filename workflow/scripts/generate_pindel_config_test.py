#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Martin Rippin"
__copyright__ = "Copyright 2021, Martin Rippin"
__email__ = "martin.rippin@igp.uu.se"
__license__ = "GPL-3"


def test_insertSize():
    from generate_pindel_config import insertSize 
    metrics = "../../.tests/integration/qc/picard_collect_multiple_metrics/HD832.HES45_T.insert_size_metrics"
    assert insertSize(metrics) == "200"
