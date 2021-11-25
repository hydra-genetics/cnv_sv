# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule svdb_merge:
    input:
        #cnvkit_vcf="cnv_sv/cnvkit_vcf/{sample}_{type}.vcf",
        #gatk_cnv_vcf="cnv_sv/gatk_cnv_vcf/{sample}_{type}.vcf",
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
        "cnv_sv/svdb_merge/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/svdb_merge/{sample}_{type}.benchmark.tsv",
            config.get("svdb_merge", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("svdb_merge", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("svdb_merge", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("svdb_merge", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("svdb_merge", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("svdb_merge", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("svdb_merge", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("svdb_merge", {}).get("container", config["default_container"])
    conda:
        "../envs/svdb_merge.yaml"
    message:
        "{rule}: Merges vcf files from different cnv callers into cnv_sv/svdb_merge/{wildcards.sample}_{wildcards.type}.merged.vcf"
    shell:
        "(svdb --merge --vcf {input.vcfs} > {output.vcf}) 2> {log}"
