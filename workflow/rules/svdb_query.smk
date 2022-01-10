# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule svdb_query:
    input:
        vcf="cnv_sv/svdb_merge/{sample}_{type}.merged.vcf",
        svdb_vcf=config.get("reference", {}).get("svdb_vcf", ""),
    output:
        vcf=temp("cnv_sv/svdb_query/{sample}_{type}.svdb_query.vcf"),
    params:
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0][:-6],
        extra=config.get("svdb_query", {}).get("extra", ""),
    log:
        "cnv_sv/svdb_query/{sample}_{type}.log",
    benchmark:
        repeat(
            "cnv_sv/svdb_query/{sample}_{type}.benchmark.tsv",
            config.get("svdb_query", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("svdb_query", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("svdb_query", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("svdb_query", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("svdb_query", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("svdb_query", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("svdb_query", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("svdb_query", {}).get("container", config["default_container"])
    conda:
        "../envs/svdb_query.yaml"
    message:
        "{rule}: Use svdb database to filter cnvs in cnv_sv/svdb_query/{wildcards.sample}_{wildcards.type}.normal_filterered.vcf"
    shell:
        "(svdb --query --query_vcf {input.vcf} --db {input.svdb_vcf} --prefix {params.prefix} {params.extra}) &> {log}"
