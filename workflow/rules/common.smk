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
if config["pacbio_alignment"] or config["ont_alignment"]:
    units = (
        pandas.read_table(config["units"], dtype=str)
        .set_index(["sample", "type", "processing_unit", "barcode"], drop=False)
        .sort_index()
    )
else:
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
            f" --mapping-bias-file={mapping_bias_file}" if mapping_bias_file is not None else "",
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


def generate_automap_id(wildcards, input):
    """
    Extracts the read group line from a bam file using pysam.

    Args:
        input: Path to the bam file.

    Returns:
        str: The read group line (e.g. @RG\\tID:group1\\tLB:library1\\tPU:unit1).
            If no read group is found, an empty string is returned.
    """

    with pysam.AlignmentFile(input.query, "rb", check_sq=False) as bam:
        # Get the header dictionary
        header = bam.header
        # Check if Read Groups are present
        if "RG" in header:
            # Access the first read group (assuming single RG in the bam)
            read_group = header["RG"][0]
            rg_line = "-R '@RG\\t" + "\\t".join(f"{key}:{val}" for key, val in read_group.items()) + "'"
            return rg_line
        else:
            return ""


def compile_output_list(wildcards):
    output_files=[] 
    #print("Check 1", config["pacbio_alignment"])
    if config["pacbio_alignment"]:
        #print ("Check 2")

        # Outputfiles which only needs proband
        files = {
            "upd": ["upd_regions.bed", "upd_sites.bed"],
        }
        output_files += [
            "cnv_sv/%s/%s_%s.%s" % (prefix, sample, type, suffix)
            for prefix in files.keys()
            for sample in samples[samples.trio_member == "proband"].index
            for type in units.loc[sample, 'type']
            for suffix in files[prefix]
        ]
        files = {
            "automap": ["HomRegions.tsv", "HomRegions.pdf"],
        }        
        output_files += [
            "cnv_sv/%s/%s_%s/%s_%s.%s" % (prefix, sample, type, sample2, type2, suffix)
            for prefix in files.keys()
            for sample in samples[samples.trio_member == "proband"].index
            for type in units.loc[sample, 'type']
            for sample2 in samples[samples.trio_member == "proband"].index
            for type2 in units.loc[sample, 'type']
            for suffix in files[prefix]
        ]
        # All output files
        files = {
            #"expansionhunter": ["vcf"],
            "pbsv_discover": ["svsig.gz"],
            "pbsv_call": ["vcf"],
        }
        output_files += [
            "cnv_sv/%s/%s_%s.%s" % (prefix, sample, type, suffix)
            for prefix in files.keys()
            for sample in samples.index  # Include all trio members
            for type in units.loc[sample, 'type']
            for suffix in files[prefix]
        ]
        #print("Outflies5", output_files)
    else:
        files = {
            "cnv_sv/cnvkit_call": ["pathology.loh.cns"],
            "cnv_sv/cnvkit_diagram": ["pdf"],
            "cnv_sv/cnvkit_scatter": ["png"],
            "cnv_sv/cnvkit_vcf": ["pathology.vcf"],
            "cnv_sv/cnvpytor": ["vcf"],
            "cnv_sv/expansionhunter": ["vcf"],
            "cnv_sv/gatk_vcf": ["pathology.vcf"],
            "cnv_sv/svdb_merge": ["no_tc.merged.vcf", "pathology.merged.vcf"],
            "cnv_sv/svdb_query": ["no_tc.svdb_query.vcf", "pathology.svdb_query.vcf"],
            "cnv_sv/exomedepth_call": ["txt", "RData"],
            "cnv_sv/pindel_vcf": ["no_tc.vcf"],
            "cnv_sv/tiddit": ["vcf"],
        }
        output_files = [
            "%s/%s_%s.%s" % (prefix, sample, unit_type, suffix)
            for prefix in files.keys()
            for sample in get_samples(samples[pd.isnull(samples["trioid"])])
            for unit_type in get_unit_types(units, sample)
            for suffix in files[prefix]
        ]
        output_files += [
            "cnv_sv/reviewer/%s_%s/" % (sample, unit_type)
            for sample in get_samples(samples[pd.isnull(samples["trioid"])])
            for unit_type in get_unit_types(units, sample)
        ]
        output_files += [
            "cnv_sv/automap/%s_%s/%s_%s.HomRegions.tsv" % (sample, unit_type, sample, unit_type)
            for sample in get_samples(samples[pd.isnull(samples["trioid"])])
            for unit_type in get_unit_types(units, sample)
        ]
        output_files.append(
            [
                "cnv_sv/manta_run_workflow_tn/%s/results/variants/somaticSV.vcf.gz" % (sample)
                for sample in get_samples(samples[pd.isnull(samples["trioid"])])
            ]
        )
        output_files.append(
            [
                "cnv_sv/manta_run_workflow_t/%s/results/variants/tumorSV.vcf.gz" % (sample)
                for sample in get_samples(samples[pd.isnull(samples["trioid"])])
            ]
        )
        output_files.append(
            [
                "cnv_sv/manta_run_workflow_n/%s/results/variants/candidateSV.vcf.gz" % (sample)
                for sample in get_samples(samples[pd.isnull(samples["trioid"])])
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
    print("Check4", output_files)
    return output_files


    