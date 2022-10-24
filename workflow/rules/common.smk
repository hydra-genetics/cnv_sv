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

min_version("7.8.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    type="N|T|R",

def get_peddy_sex(wildcards, peddy_sex_check):
    sample = '{}_{}'.format(wildcards.sample, wildcards.type)
    sex_df = pd.read_table(peddy_sex_check, sep=',').set_index("sample_id", drop=False)

    sample_sex = sex_df.at[sample, 'predicted_sex']

    return sample_sex


def get_locus_str(loci):
    with open(loci, 'r') as catfile:
        loc_str = catfile.readline().rstrip()
    return loc_str


def compile_output_list(wildcards):
    files = {
        "cnv_sv/cnvkit_call": ["loh.cns"],
        "cnv_sv/cnvkit_diagram": ["pdf"],
        "cnv_sv/cnvkit_scatter": ["png"],
        "cnv_sv/cnvkit_vcf": ["vcf"],
        "cnv_sv/expansionhunter": ["vcf"],
        "cnv_sv/gatk_vcf": ["vcf"],
        "cnv_sv/svdb_merge": ["merged.vcf"],
        "cnv_sv/svdb_query": ["svdb_query.vcf"],
        "cnv_sv/exomedepth_call": ["SV.txt"],
        "cnv_sv/pindel_vcf": ["vcf"],
        "cnv_sv/tiddit": ["vcf"],
    }
    output_files = [
        "%s/%s_%s.%s" % (prefix, sample, unit_type, suffix)
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for suffix in files[prefix]
    ]
    output_files += [
        "cnv_sv/expansionhunter/reviewer/%s_%s/" % (sample, unit_type)
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]
    output_files.append(
        ["cnv_sv/manta_run_workflow_tn/%s/results/variants/somaticSV.vcf.gz" % (sample) for sample in get_samples(samples)]
    )
    output_files.append(
        ["cnv_sv/manta_run_workflow_t/%s/results/variants/tumorSV.vcf.gz" % (sample) for sample in get_samples(samples)]
    )

    output_files.append(
        ["cnv_sv/manta_run_workflow_n/%s/results/variants/candidateSV.vcf.gz" % (sample) for sample in get_samples(samples)]
    )

    return output_files
