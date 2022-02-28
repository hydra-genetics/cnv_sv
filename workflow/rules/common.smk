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
    files = {
        "cnv_sv/cnvkit_call": [
            "loh.cns",
        ],
        "cnv_sv/cnvkit_diagram": [
            "pdf",
        ],
        "cnv_sv/cnvkit_scatter": [
            "png",
        ],
        "cnv_sv/cnvkit_vcf": [
            "vcf",
        ],
        "cnv_sv/gatk_cnv_vcf": [
            "vcf",
        ],
        "cnv_sv/svdb_merge": [
            "merged.vcf",
        ],
        "cnv_sv/svdb_query": [
            "svdb_query.vcf",
        ],
    }
    output_files = [
        "%s/%s_%s.%s" % (prefix, sample, unit_type, suffix)
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for suffix in files[prefix]
    ]
    output_files.append(
        ["cnv_sv/pindel/%s.vcf" % (sample) for sample in get_samples(samples)]
    )
    output_files.append(
        ["cnv_sv/manta_run_workflow_tn/%s/results/variants/somaticSV.vcf.gz" % (sample) for sample in get_samples(samples)]
    )
    output_files.append(
        ["cnv_sv/manta_run_workflow_t/%s/results/variants/tumorSV.vcf.gz" % (sample) for sample in get_samples(samples)]
    )
    return output_files
