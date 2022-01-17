# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule gatk_cnv_collect_allelic_counts:
    input:
        bam="alignment/merge_bam/{sample}_{type}.bam",
        bai="alignment/merge_bam/{sample}_{type}.bam.bai",
        interval=config.get("gatk_cnv_collect_allelic_counts", {}).get("SNP_interval", ""),
        ref=config["reference"]["fasta"],
    output:
        temp("cnv_sv/gatk_cnv_collect_allelic_counts/{sample}_{type}.clean.allelicCounts.tsv"),
    params:
        extra=config.get("gatk_cnv_collect_allelic_counts", {}).get("extra", ""),
    log:
        "cnv_sv/gatk_cnv_collect_allelic_counts/{sample}_{type}.clean.allelicCounts.tsv.log",
    benchmark:
        repeat(
            "cnv_sv/gatk_cnv_collect_allelic_counts/{sample}_{type}.clean.allelicCounts.tsv.benchmark.tsv",
            config.get("gatk_cnv_collect_allelic_counts", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_cnv_collect_allelic_counts", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("gatk_cnv_collect_allelic_counts", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_cnv_collect_allelic_counts", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("gatk_cnv_collect_allelic_counts", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_cnv_collect_allelic_counts", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("gatk_cnv_collect_allelic_counts", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("gatk_cnv_collect_allelic_counts", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk_cnv_collect_allelic_counts.yaml"
    message:
        "{rule}: Use gatk_cnv to obtain cnv_sv/gatk_cnv_collect_allelic_counts/{wildcards.sample}_{wildcards.type}.clean.allelicCounts.tsv"
    shell:
        "(gatk --java-options '-Xmx4g' CollectAllelicCounts "
        "-L {input.interval} "
        "-I {input.bam} "
        "-R {input.ref} "
        "-O {output}"
        "{params.extra}) &> {log}"
