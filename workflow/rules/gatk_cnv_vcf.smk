# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule gatk_cnv_vcf:
    input:
        segment="cnv_sv/gatk_cnv_model_segments/{sample}_{type}.clean.modelFinal.seg",
    output:
        vcf=temp("cnv_sv/gatk_cnv_vcf/{sample}_{type}.vcf"),
    params:
        sample_id="{sample}_{type}",
        hom_del_limit=config.get("gatk_cnv_vcf", {}).get("hom_del_limit", 0.5),
        het_del_limit=config.get("gatk_cnv_vcf", {}).get("het_del_limit", 1.5),
        dup_limit=config.get("gatk_cnv_vcf", {}).get("dup_limit", 2.5),
        TC=lambda wildcards: get_sample(samples, wildcards)["tumor_content"],
    log:
        "cnv_sv/gatk_cnv_vcf/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/gatk_cnv_vcf/{sample}_{type}.vcf.benchmark.tsv",
            config.get("gatk_cnv_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_cnv_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("gatk_cnv_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_cnv_vcf", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("gatk_cnv_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_cnv_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_cnv_vcf", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("gatk_cnv_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk_cnv_vcf.yaml"
    message:
        "{rule}: Export gatk cnv segments into vcf in cnv_sv/gatk_cnv_vcf/{wildcards.sample}_{wildcards.type}.vcf"
    script:
        "../scripts/gatk_cnv_vcf.py"
