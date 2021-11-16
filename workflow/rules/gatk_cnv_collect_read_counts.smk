# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule gatk_cnv_collect_read_counts:
    input:
        bam="alignment/merge_bam/{sample}_{type}.bam",
        bai="alignment/merge_bam/{sample}_{type}.bam.bai",
        interval=config["gatk_cnv_collect_read_counts"]["interval"],
    output:
        temp("cnv_sv/gatk_cnv_collect_read_counts/{sample}_{type}.counts.hdf5"),
    params:
        mergingRule="OVERLAPPING_ONLY",
        extra=config.get("gatk_cnv_collect_read_counts", {}).get("extra", ""),
    log:
        "cnv_sv/gatk_cnv_collect_read_counts/{sample}_{type}.counts.hdf5.log",
    benchmark:
        repeat(
            "cnv_sv/gatk_cnv_collect_read_counts/{sample}_{type}.counts.hdf5.benchmark.tsv",
            config.get("gatk_cnv_collect_read_counts", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_cnv_collect_read_counts", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("gatk_cnv_collect_read_counts", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_cnv_collect_read_counts", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("gatk_cnv_collect_read_counts", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_cnv_collect_read_counts", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_cnv_collect_read_counts", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("gatk_cnv_collect_read_counts", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk_cnv_collect_read_counts.yaml"
    message:
        "{rule}: Use gatk_cnv to obtain cnv_sv/gatk_cnv_collect_read_counts/{wildcards.sample}_{wildcards.type}.counts.hdf5"
    shell:
        "(gatk --java-options '-Xmx4g' CollectReadCounts "
        "-I {input.bam} "
        "-L {input.interval} "
        "--interval-merging-rule {params.mergingRule} "
        "{params.extra} "
        "-O {output}) &> {log}"
