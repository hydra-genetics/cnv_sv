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


def get_purecn_inputs(wildcards: snakemake.io.Wildcards):
    inputs = {k: v for k, v in config.get("purecn", {}).items() if k in ["normaldb", "mapping_bias_file", "snp_blacklist"]}
    segmentation_method = config.get("purecn", {}).get("segmentation_method", "")
    if segmentation_method == "internal":
        inputs.update(
            {
                "tumor": f"cnv_sv/purecn_coverage/{wildcards.sample}_T_coverage_loess.txt.gz",
                "intervals": config.get("purecn", {}).get("intervals"),
                "normaldb": config.get("purecn", {}).get("normaldb"),
            }
        )
    elif segmentation_method == "GATK4":
        inputs.update(
            {
                "tumor": f"cnv_sv/gatk_collect_read_counts/{wildcards.sample}_T.counts.hdf5",
                "seg_file": f"cnv_sv/gatk_model_segments/{wildcards.sample}_T.clean.modelFinal.seg",
                "log_ratio_file": f"cnv_sv/gatk_denoise_read_counts/{wildcards.sample}_T.clean.denoisedCR.tsv",
            }
        )
    elif segmentation_method == "cnvkit":
        inputs.update(
            {
                "tumor": f"cnv_sv/cnvkit_batch/{wildcards.sample}/{wildcards.sample}_T.cnr",
                "seg_file": f"cnv_sv/cnvkit_export_seg/{wildcards.sample}_T.seg",
            }
        )

    return inputs


def get_purecn_extra(wildcards: snakemake.io.Wildcards, input: snakemake.io.InputFiles, threads: int):
    log_ratio_file = input.get("log_ratio_file")
    seg_file = input.get("seg_file")
    intervals = input.get("intervals")
    normaldb = input.get("normaldb")
    mapping_bias_file = input.get("mapping_bias_file")
    snp_blacklist = input.get("snp_blacklist")

    extra = "".join(
        [
            config.get("purecn", {}).get("extra", ""),
            f" --log-ratio-file={log_ratio_file}" if log_ratio_file is not None else "",
            f" --seg-file={seg_file}" if seg_file is not None else "",
            f" --intervals={intervals}" if intervals is not None else "",
            f" --mapping-bias-file={mapping_bias_file}" if mapping_bias_file is not None else "",
            f" --normaldb={normaldb}" if normaldb is not None else "",
            f" --snp-blacklist={snp_blacklist}" if snp_blacklist is not None else "",
            f" --parallel --cores={threads}" if threads > 1 else "",
        ]
    )

    return extra


def compile_output_list(wildcards):
    output_files = {
        path
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for path in [
            f"cnv_sv/cnvkit_call/{sample}_{unit_type}.loh.cns",
            f"cnv_sv/cnvkit_diagram/{sample}_{unit_type}.pdf",
            f"cnv_sv/cnvkit_scatter/{sample}_{unit_type}.png",
            f"cnv_sv/cnvkit_vcf/{sample}_{unit_type}.vcf",
            f"cnv_sv/cnvkit_seg/{sample}_{unit_type}.seg",
            f"cnv_sv/gatk_vcf/{sample}_{unit_type}.vcf",
            f"cnv_sv/svdb_merge/{sample}_{unit_type}.merged.vcf",
            f"cnv_sv/svdb_query/{sample}_{unit_type}.svdb_query.vcf",
            f"cnv_sv/exomedepth_call/{sample}_{unit_type}.SV.txt",
            f"cnv_sv/pindel_vcf/{sample}_{unit_type}.vcf",
            f"cnv_sv/manta_run_workflow_tn/{sample}/results/variants/somaticSV.vcf.gz",
            f"cnv_sv/purecn_coverage/{sample}_{unit_type}_coverage_loess.txt.gz",
            f"cnv_sv/purecn/{sample}_T.csv",
        ]
    }

    return [*output_files]
