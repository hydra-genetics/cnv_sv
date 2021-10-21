# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule GATK_cnv_collectAllelicCounts:
    input:
        bam="alignment/bwa_mem/{sample}_{type}.bam",
        bai="alignment/bwa_mem/{sample}_{type}.bam.bai",
        interval=config["GATK_cnv_collectAllelicCounts"]["SNP_interval"],
        ref=config["reference"]["fasta"],
    output:
        temp("cnv/GATK_cnv_collectAllelicCounts/{sample}_{type}.clean.allelicCounts.tsv"),
    params:
        mergingRule="OVERLAPPING_ONLY",
        extra=config.get("GATK_cnv_collectAllelicCounts", {}).get("extra", ""),
    log:
        "cnv/GATK_cnv_collectAllelicCounts/{sample}_{type}.clean.allelicCounts.tsv.log",
    benchmark:
        repeat(
            "cnv/GATK_cnv_collectAllelicCounts/{sample}_{type}.clean.allelicCounts.tsv.benchmark.tsv",
            config.get("GATK_cnv_collectAllelicCounts", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("GATK_cnv_collectAllelicCounts", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("GATK_cnv_collectAllelicCounts", {}).get("container", config["default_container"])
    conda:
        "../envs/GATK_cnv_collectAllelicCounts.yaml"
    message:
        "{rule}: Use GATK_cnv to obtain cnv/GATK_cnv_collectAllelicCounts/{wildcards.sample}_{wildcards.type}.clean.allelicCounts.tsv"
    shell:
        "(gatk --java-options '-Xmx4g' CollectAllelicCounts "
        "-L {input.interval} "
        "-I {input.bam} "
        "-R {input.ref} "
        "-O {output}"
        "{params.extra}) &> {log}"
