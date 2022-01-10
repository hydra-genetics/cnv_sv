# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

min_version("6.8.0")

### Set and validate config file


configfile: "config.yaml"


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pandas.read_table(config["units"], dtype=str).set_index(["sample", "type", "run", "lane"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    unit="N|T|R",


def compile_output_list(wildcards):
    output_files = [
        "cnv_sv/cnvkit_call_loh/%s_%s.loh.cns" % (sample, t)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
    ]
    output_files.append(
        [
            "cnv_sv/gatk_cnv_call_copy_ratio_segments/%s_%s.clean.calledCNVs.seg" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        ["cnv_sv/gatk_cnv_vcf/%s_%s.vcf" % (sample, t) for sample in get_samples(samples) for t in get_unit_types(units, sample)]
    )
    output_files.append(
        ["cnv_sv/cnvkit_vcf/%s_%s.vcf" % (sample, t) for sample in get_samples(samples) for t in get_unit_types(units, sample)]
    )
    output_files.append(
        [
            "cnv_sv/svdb_merge/%s_%s.merged.vcf" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "cnv_sv/svdb_query/%s_%s.svdb_query.vcf" % (sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    )
    return output_files
