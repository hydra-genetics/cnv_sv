# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule GATK_cnv_collectReadCounts:
    input:
        bam="alignment/bwa_mem/{sample}_{type}.bam",
        bai="alignment/bwa_mem/{sample}_{type}.bam.bai",
        interval=config["GATK_cnv_collectReadCounts"]["interval"],
    output:
        "cnv/GATK_cnv_collectReadCounts/{sample}_{type}.counts.hdf5",
    params:
        mergingRule="OVERLAPPING_ONLY",
        extra=config.get("GATK_cnv_collectReadCounts", {}).get("extra", ""),
    log:
        "cnv/GATK_cnv_collectReadCounts/{sample}_{type}.counts.hdf5.log",
    benchmark:
        repeat(
            "cnv/GATK_cnv_collectReadCounts/{sample}_{type}.counts.hdf5.benchmark.tsv",
            config.get("GATK_cnv_collectReadCounts", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("GATK_cnv_collectReadCounts", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("GATK_cnv_collectReadCounts", {}).get("container", config["default_container"])
    conda:
        "../envs/GATK_cnv_collectReadCounts.yaml"
    message:
        "{rule}: Use GATK_cnv to obtain cnv/GATK_cnv_collectReadCounts/{wildcards.sample}_{wildcards.type}.counts.hdf5"
    shell:
        "(gatk --java-options '-Xmx4g' CollectReadCounts "
        "-I {input.bam} "
        "-L {input.interval} "
        "--interval-merging-rule {params.mergingRule} "
        "{params.extra} "
        "-O {output}) &> {log}"
