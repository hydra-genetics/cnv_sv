# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule svdb_merge:
    input:
        vcfs=expand(
            "cnv_sv/{cnv_caller}_vcf/{{sample}}_{{type}}.vcf",
            cnv_caller=config.get("svdb_merge", {}).get("cnv_callers", []),
        ),
    output:
        vcf=temp("cnv_sv/svdb_merge/{sample}_{type}.merged.vcf"),
    params:
        overlap=config.get("svdb_merge", {}).get("overlap", 0.6),
        extra=config.get("svdb_merge", {}).get("extra", ""),
    log:
        "cnv_sv/svdb_merge/{sample}_{type}.merged.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/svdb_merge/{sample}_{type}.merged.benchmark.tsv",
            config.get("svdb_merge", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("svdb_merge", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("svdb_merge", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("svdb_merge", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("svdb_merge", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("svdb_merge", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("svdb_merge", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("svdb_merge", {}).get("container", config["default_container"])
    conda:
        "../envs/svdb.yaml"
    message:
        "{rule}: Merges vcf files from different cnv callers into {output.vcf}"
    shell:
        "(svdb --merge --vcf {input.vcfs} > {output.vcf}) 2> {log}"


rule svdb_query:
    input:
        vcf="cnv_sv/svdb_merge/{sample}_{type}.merged.vcf",
        svdb_vcf=config.get("svdb_query", {}).get("svdb_vcf", ""),
    output:
        vcf=temp("cnv_sv/svdb_query/{sample}_{type}.svdb_query.vcf"),
    params:
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0][:-6],
        extra=config.get("svdb_query", {}).get("extra", ""),
    log:
        "cnv_sv/svdb_query/{sample}_{type}.svdb_query.log",
    benchmark:
        repeat(
            "cnv_sv/svdb_query/{sample}_{type}.svdb_query.benchmark.tsv",
            config.get("svdb_query", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("svdb_query", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("svdb_query", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("svdb_query", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("svdb_query", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("svdb_query", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("svdb_query", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("svdb_query", {}).get("container", config["default_container"])
    conda:
        "../envs/svdb.yaml"
    message:
        "{rule}: Use svdb database to filter cnvs into {output.vcf}"
    shell:
        "(svdb --query --query_vcf {input.vcf} --db {input.svdb_vcf} --prefix {params.prefix} {params.extra}) &> {log}"
