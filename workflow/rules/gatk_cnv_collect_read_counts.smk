# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule gatk_cnv_collectReadCounts:
    input:
        bam="alignment/bwa_mem/{sample}_{type}.bam",
        bai="alignment/bwa_mem/{sample}_{type}.bam.bai",
        interval=config["gatk_cnv_collect_read_counts"]["interval"],
    output:
        temp("cnv/gatk_cnv_collect_read_counts/{sample}_{type}.counts.hdf5"),
    params:
        mergingRule="OVERLAPPING_ONLY",
        extra=config.get("gatk_cnv_collect_read_counts", {}).get("extra", ""),
    log:
        "cnv/gatk_cnv_collect_read_counts/{sample}_{type}.counts.hdf5.log",
    benchmark:
        repeat(
            "cnv/gatk_cnv_collect_read_counts/{sample}_{type}.counts.hdf5.benchmark.tsv",
            config.get("gatk_cnv_collect_read_counts", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_cnv_collect_read_counts", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("gatk_cnv_collect_read_counts", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk_cnv_collect_read_counts.yaml"
    message:
        "{rule}: Use gatk_cnv to obtain cnv/gatk_cnv_collect_read_counts/{wildcards.sample}_{wildcards.type}.counts.hdf5"
    shell:
        "(gatk --java-options '-Xmx4g' CollectReadCounts "
        "-I {input.bam} "
        "-L {input.interval} "
        "--interval-merging-rule {params.mergingRule} "
        "{params.extra} "
        "-O {output}) &> {log}"
