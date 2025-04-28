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

min_version("7.8.3")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype={"sample": str}).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file
units = pandas.read_table(config["units"], dtype=str)

if units.platform.iloc[0] in ["PACBIO", "ONT"]:
    units = units.set_index(["sample", "type", "processing_unit", "barcode"], drop=False).sort_index()
else:  # assume that the platform Illumina data with a lane and flowcell columns
    units = units.set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    type="N|T|R",
    file="^cnv_sv/.+",


def get_longread_bam(wildcards):
    aligner = config.get("aligner", "minimap2")
    alignment_path = f"alignment/{aligner}_align/{wildcards.sample}_{wildcards.type}.bam"
    index_path = f"alignment/{aligner}_align/{wildcards.sample}_{wildcards.type}.bam.bai"
    return alignment_path, index_path


def get_input_bam(wildcards, default_path="alignment/samtools_merge_bam"):
    """
    Get path to input bam files.

    This function checks if 'haplotagged_bam', 'aligner' and 'haplotagging' are in the config,
    interprets their combination to compile paths to input bam files.

    The 'aligner' entry should contain name of the aligner used. For
    example, if the aligner used is minimap2, the entry should be

        "aligner": "minimap2",

    The 'haplotagged_bam' entry should be set to true in the config if the bam input is already haplotagged.
    The entry should look like this:

        "haplotagged_bam": true,

    If the input bam files are not haplotagged, the entry should be set to false:
        "haplotagged_bam": false,

    The 'haplotagging_tool' entry should contain name of the haplotagging tool used.
    For example, if the tool used is whatshap, the entry should be
        "haplotagging": "whatshap",

    If neither 'haplotagged_bam' nor 'aligner' are in the config, it defaults to
    'alignment/samtools_merge_bam'


    The function returns path to the input bam file and its index file.

    Arguments:
    wildcards: snakemake.io.Wildcards
        The wildcards object containing the sample name and type.

    Returns:
    tuple: (alignment_path, index_path)
        The path to the input bam file and its index file.
    """

    if config.get("haplotagged_bam") is True and config.get("aligner") is None:
        # filter lines with status 'haplotagged' from units.tsv
        # use the string from column 'bam' as input path
        unit = units[
            (units["sample"] == wildcards.sample) & (units["type"] == wildcards.type) & (units["status"] == "haplotagged")
        ]
        alignment_path = unit["bam"].iloc[0]
        index_path = f"{alignment_path}.bai"

    elif config.get("haplotagged_bam") is None and config.get("aligner") is not None:
        # if haplotagged_bam entry is not in the config & aligner is in the config
        # use get_longread_bam to get bam paths
        path_index = get_longread_bam(wildcards)
        alignment_path = path_index[0]
        index_path = path_index[1]

    elif config.get("haplotagged_bam") is False and config.get("haplotagging") is not None:
        # if config contains haplotagging_tool and haplotagged_bam: false
        # use this tool to compile input bam path
        tool = config.get("haplotagging")
        alignment_path = f"annotation/{tool}_haplotag/{wildcards.sample}_{wildcards.type}.haplotagged.bam"
        index_path = f"annotation/{tool}_haplotag/{wildcards.sample}_{wildcards.type}.haplotagged.bam.bai"

    else:
        # if neither input_bam nor aligner entries are in the config
        # use default bam path to compile output
        alignment_path = f"{default_path}/{wildcards.sample}_{wildcards.type}.bam"
        index_path = f"{default_path}/{wildcards.sample}_{wildcards.type}.bam.bai"
    return alignment_path, index_path


def get_karyotype(wildcards):
    """
    Translate sex to karyotype for trgt. If sex is unknown
    or sex is not in the samples.tsv file then 'XX' is returned
    as default
    """

    sex = samples.loc[wildcards.sample].sex
    if sex == "male":
        karyotype = "XY"
    else:
        karyotype = "XX"

    return karyotype


def get_expected_cn(wildcards):
    """
    Find the expected copy number BED file path for sawfish discover.
    These bed are typically used to specify ploidy in the non-PAR regions of the sex chromosomes.
    """

    sex = samples.loc[wildcards.sample].sex
    if sex == "male":
        expected_cn = config.get("sawfish_discover", {}).get("expected_cn", {}).get("male", "")
        sawfish_param = f"--expected-cn {expected_cn}"
    elif sex == "female":
        expected_cn = config.get("sawfish_discover", {}).get("expected_cn", {}).get("female", "")
        sawfish_param = f"--expected-cn {expected_cn}"
    else:  # when no sex in samples.tsv treat all regions as diploid
        sawfish_param = ""

    return sawfish_param


def get_tc(wildcards):
    """
    Get the tumor cell content of a sample. If the tumor content
    cannot be identified, return an empty string.
    """
    tc_method = wildcards.tc_method
    if tc_method == "pathology":
        sample = get_sample(samples, wildcards)
        if not "tumor_content" in sample:
            return ""
        return get_sample(samples, wildcards)["tumor_content"]
    else:
        tc_file = f"cnv_sv/{tc_method}_purity_file/{wildcards.sample}_{wildcards.type}.purity.txt"
        if not os.path.exists(tc_file):
            return ""
        else:
            with open(tc_file, "r") as f:
                tc = f.read()
            return tc


def get_tc_file(wildcards):
    tc_method = wildcards.tc_method
    if tc_method == "pathology":
        return config.get("samples")
    else:
        return f"cnv_sv/{tc_method}_purity_file/{wildcards.sample}_{wildcards.type}.purity.txt"


def get_purecn_inputs(wildcards: snakemake.io.Wildcards):
    inputs = {k: v for k, v in config.get("purecn", {}).items() if k in ["normaldb", "mapping_bias_file", "snp_blacklist"]}
    segmentation_method = config.get("purecn", {}).get("segmentation_method", "")
    if segmentation_method == "internal":
        inputs.update(
            {
                "tumor": f"cnv_sv/purecn_coverage/{wildcards.sample}_{wildcards.type}_coverage_loess.txt.gz",
                "intervals": config.get("purecn", {}).get("intervals"),
                "normaldb": config.get("purecn", {}).get("normaldb"),
            }
        )
    elif segmentation_method == "GATK4":
        inputs.update(
            {
                "tumor": f"cnv_sv/gatk_collect_read_counts/{wildcards.sample}_{wildcards.type}.counts.hdf5",
                "seg_file": f"cnv_sv/gatk_model_segments/{wildcards.sample}_{wildcards.type}.clean.modelFinal.seg",
                "log_ratio_file": f"cnv_sv/gatk_denoise_read_counts/{wildcards.sample}_{wildcards.type}.clean.denoisedCR.tsv",
            }
        )
    elif segmentation_method == "cnvkit":
        inputs.update(
            {
                "tumor": f"cnv_sv/cnvkit_batch/{wildcards.sample}/{wildcards.sample}_{wildcards.type}.cnr",
                "seg_file": f"cnv_sv/cnvkit_export_seg/{wildcards.sample}_{wildcards.type}.seg",
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
            (f" --mapping-bias-file={mapping_bias_file}" if mapping_bias_file is not None else ""),
            f" --normaldb={normaldb}" if normaldb is not None else "",
            f" --snp-blacklist={snp_blacklist}" if snp_blacklist is not None else "",
            f" --parallel --cores={threads}" if threads > 1 else "",
        ]
    )

    return extra


def get_peddy_sex(wildcards, peddy_sex_check):
    sample = "{}_{}".format(wildcards.sample, wildcards.type)
    sex_df = pd.read_table(peddy_sex_check, sep=",").set_index("sample_id", drop=False)

    sample_sex = sex_df.at[sample, "predicted_sex"]

    return sample_sex


def get_exomedepth_ref(wildcards):
    sex = get_peddy_sex(wildcards, checkpoints.exomedepth_sex.get().output[0])

    if sex == "male":
        ref = config.get("exomedepth_call", {}).get("male_reference", "")
    else:  # use female ref in the case of female or NA
        ref = config.get("exomedepth_call", {}).get("female_reference", "")

    return ref


def get_locus_str(loci):
    with open(loci, "r") as catfile:
        loc_str = catfile.readline().rstrip()
    return loc_str


def get_vcfs_for_svdb_merge(wildcards, add_suffix=False):
    vcf_dict = {}
    for v in config.get("svdb_merge", {}).get("tc_method"):
        tc_method = v["name"]
        callers = v["cnv_caller"]
        for caller in callers:
            if add_suffix:
                caller_suffix = f":{caller}"
            else:
                caller_suffix = ""
            if tc_method in vcf_dict:
                vcf_dict[tc_method].append(
                    f"cnv_sv/{caller}_vcf/{wildcards.sample}_{wildcards.type}.{tc_method}.vcf{caller_suffix}"
                )
            else:
                vcf_dict[tc_method] = [f"cnv_sv/{caller}_vcf/{wildcards.sample}_{wildcards.type}.{tc_method}.vcf{caller_suffix}"]
    return vcf_dict[wildcards.tc_method]


def get_priority(wildcards):
    priority_dict = {}
    for v in config.get("svdb_merge", {}).get("tc_method"):
        tc_method = v["name"]
        priority_dict[tc_method] = v["priority"]

    return priority_dict[wildcards.tc_method]


def get_parent_samples(wildcards, trio_member):
    proband_sample = samples[samples.index == wildcards.sample]
    trio_id = proband_sample.at[wildcards.sample, "trioid"]

    parent_sample = samples[(samples.trio_member == trio_member) & (samples.trioid == trio_id)].index[0]

    parent_sample_id = f"{parent_sample}_{wildcards.type}"

    return parent_sample_id


def get_trgt_loci(wildcards):
    trgt_bed = config.get("trgt_genotype", {}).get("bed", "")
    rep_ids = []
    with open(trgt_bed, "r") as infile:
        for line in infile:
            cols = line.split("\t")
            rep_id = cols[3].split(";")[0]
            rep_ids.append(rep_id.split("=")[1])
    return rep_ids


def get_tr_bed(wildcards):
    tr_bed = config.get("sniffles2_call", {}).get("tandem_repeats", "")

    if tr_bed != "":
        tr_bed = f"--tandem-repeats {tr_bed}"

    return tr_bed


def compile_output_list(wildcards):
    platform = units.platform.iloc[0]
    files = {
        "cnv_sv/trgt_genotype": ["vcf.gz"],
        "cnv_sv/sniffles2_call": ["vcf.gz", "vcf.gz.tbi", "snf"],
        "cnv_sv/sawfish_joint_call": ["vcf.gz"],
    }
    output_files = [
        f"{prefix}/{sample}_{unit_type}.{suffix}"
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for platform in units.loc[(sample,)].platform
        if platform in ["ONT", "PACBIO"]
        for suffix in files[prefix]
    ]

    files = {
        "cnv_sv/trgt_plot": [config.get("trgt_plot", {}).get("image", "svg")],
    }
    output_files += [
        f"{prefix}/{sample}_{unit_type}_{locus}.{suffix}"
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for platform in units.loc[(sample,)].platform
        if platform in ["ONT", "PACBIO"]
        for locus in get_trgt_loci(wildcards)
        for suffix in files[prefix]
    ]

    files = {
        "cnv_sv/paraphase": ["bam"],
        "cnv_sv/paraphase": [".bam.bai"],
        "cnv_sv/paraphase": ["json"],
    }
    output_files += [
        f"{prefix}/paraphase_{sample}_{unit_type}/{sample}_{unit_type}.paraphase.{suffix}"
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for platform in units.loc[(sample,)].platform
        if platform in ["ONT", "PACBIO"]
        for suffix in files[prefix]
    ]

    files = {
        "cnv_sv/paraphase": ["vcf.gz"],
    }
    output_files += [
        f"{prefix}/paraphase_{sample}_{unit_type}/{sample}_{unit_type}_paraphase_vcfs/{sample}_{unit_type}_{gene}.{suffix}"
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for platform in units.loc[(sample,)].platform
        if platform in ["ONT", "PACBIO"]
        for gene in config.get("paraphase", {}).get("genes", "")
        for suffix in files[prefix]
    ]

    files = {
        "cnv_sv/cnvkit_call": ["pathology.loh.cns"],
        "cnv_sv/cnvkit_diagram": ["pdf"],
        "cnv_sv/cnvkit_scatter": ["png"],
        "cnv_sv/cnvkit_vcf": ["pathology.vcf.gz"],
        "cnv_sv/cnvpytor": ["vcf.gz"],
        "cnv_sv/expansionhunter": ["vcf.gz"],
        "cnv_sv/gatk_vcf": ["pathology.vcf.gz"],
        "cnv_sv/svdb_merge": ["no_tc.merged.vcf.gz", "pathology.merged.vcf.gz"],
        "cnv_sv/svdb_query": ["no_tc.svdb_query.vcf.gz", "pathology.svdb_query.vcf.gz"],
        "cnv_sv/exomedepth_call": ["txt", "RData"],
        "cnv_sv/pindel_vcf": ["no_tc.vcf.gz"],
        "cnv_sv/tiddit": ["vcf.gz"],
        "cnv_sv/scanitd": ["vcf"],
    }
    output_files += [
        "%s/%s_%s.%s" % (prefix, sample, unit_type, suffix)
        for prefix in files.keys()
        for sample in get_samples(samples[pd.isnull(samples["trioid"])])
        for unit_type in get_unit_types(units, sample)
        for platform in units.loc[(sample,)].platform
        if platform not in ["ONT", "PACBIO"]
        for suffix in files[prefix]
    ]
    output_files += [
        "cnv_sv/reviewer/%s_%s/" % (sample, unit_type)
        for sample in get_samples(samples[pd.isnull(samples["trioid"])])
        for unit_type in get_unit_types(units, sample)
        for platform in units.loc[(sample,)].platform
        if platform not in ["ONT", "PACBIO"]
    ]
    output_files += [
        "cnv_sv/automap/%s_%s/%s_%s.HomRegions.tsv" % (sample, unit_type, sample, unit_type)
        for sample in get_samples(samples[pd.isnull(samples["trioid"])])
        for unit_type in get_unit_types(units, sample)
        for platform in units.loc[(sample,)].platform
        if platform not in ["ONT", "PACBIO"]
    ]
    output_files.append(
        [
            "cnv_sv/manta_run_workflow_tn/%s/results/variants/somaticSV.vcf.gz" % (sample)
            for sample in get_samples(samples[pd.isnull(samples["trioid"])])
            for platform in units.loc[(sample,)].platform
            if platform not in ["ONT", "PACBIO"]
        ]
    )
    output_files.append(
        [
            "cnv_sv/manta_run_workflow_t/%s/results/variants/tumorSV.vcf.gz" % (sample)
            for sample in get_samples(samples[pd.isnull(samples["trioid"])])
            for platform in units.loc[(sample,)].platform
            if platform not in ["ONT", "PACBIO"]
        ]
    )
    output_files.append(
        [
            "cnv_sv/manta_run_workflow_n/%s/results/variants/candidateSV.vcf.gz" % (sample)
            for sample in get_samples(samples[pd.isnull(samples["trioid"])])
            for platform in units.loc[(sample,)].platform
            if platform not in ["ONT", "PACBIO"]
        ]
    )
    files = {
        "upd": ["upd_regions.bed", "upd_sites.bed"],
    }
    output_files += [
        "cnv_sv/%s/%s_%s.%s" % (prefix, sample, unit_type, suffix)
        for prefix in files.keys()
        for sample in samples[samples.trio_member == "proband"].index
        for unit_type in get_unit_types(units, sample)
        for platform in units.loc[(sample,)].platform
        if platform not in ["ONT", "PACBIO"]
        for suffix in files[prefix]
    ]
    files = {
        "cnv_sv/pbsv_discover": ["svsig.gz"],
        "cnv_sv/pbsv_call": ["vcf"],
    }
    output_files += [
        f"{prefix}/{sample}_{unit_type}.{suffix}"
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for platform in units.loc[(sample,)].platform
        if platform in ["ONT", "PACBIO"]
        for suffix in files[prefix]
    ]
    files = {
        "cnv_sv/hificnv": ["vcf.gz"],
        "cnv_sv/hificnv": ["depth.bw"],
        "cnv_sv/hificnv": ["copynum.bedgraph"],
    }
    output_files += [
        f"{prefix}/{sample}_{unit_type}.{suffix}"
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for platform in units.loc[(sample,)].platform
        if platform in ["ONT", "PACBIO"]
        for suffix in files[prefix]
    ]

    # Since it is not possible to create integration test without a full dataset purecn will not be subjected to integration
    # testing and we can not guarantee that it will work
    # output_files.append(
    #   [
    #       "cnv_sv/purecn_coverage/%s_%s_coverage_loess.txt.gz" % (sample, unit_type)
    #       for sample in get_samples(samples)
    #       for unit_type in get_unit_types(units, sample)
    #   ]
    # )
    # output_files.append(["cnv_sv/purecn/%s_T.csv" % (sample) for sample in get_samples(samples)])

    # Since it is not possible to create integration test without a large dataset SMNCopyNumberCaller will not be subjected to integration
    # testing and we can not guarantee that it will work
    # output_files += [
    #     ""cnv_sv/smn_caller"" % (sample, unit_type)
    #     for sample in get_samples(samples)
    #     for unit_type in get_unit_types(units, sample)
    # ]
    # output_files += [
    #     "cnv_sv/smn_charts/smn_%s_%s.pdf" % (sample, unit_type)
    #     for sample in get_samples(samples)
    #     for unit_type in get_unit_types(units, sample)
    # ]

    # Since it is not possible to create integration test without a large dataset jumble will not be subjected to integration
    # testing and we can not guarantee that it will work
    # output_files += [
    #     "cnv_sv/jumble_vcf/%s_%s.pathology.vcf" % (sample, unit_type)
    #     for sample in get_samples(samples)
    #     for unit_type in get_unit_types(units, sample)
    # ]

    # Can't access the newest version of MELT for integration-test right now. Add later when we have docker with newest version of MELT.
    # files = {
    #    "melt": ["ALU.final_comp.vcf", "LINE1.final_comp.vcf", "SVA.final_comp.vcf", HERVK.final_comp.vcf],
    # }
    # output_files += [
    #    "cnv_sv/%s/%s_%s/%s" % (prefix, sample, unit_type, suffix)
    #    for prefix in files.keys()
    #    for sample in get_samples(samples[pd.isnull(samples["trioid"])])
    #    for unit_type in get_unit_types(units, sample)
    #    for suffix in files[prefix]
    # ]

    return output_files
