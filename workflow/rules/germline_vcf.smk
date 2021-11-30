# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule germline_vcf:
    input:
        vcf="snv_indels/ensemble_vcf/{sample}_{type}.ensembled.vep_annotated.vcf",
    output:
        vcf=temp("cnv_sv/germline_vcf/{sample}_{type}.germline.vcf"),
    params:
        filter=config.get("germline_vcf", {}).get(
            "filter",
            '--filter "DP > 50" --filter "AF >= 0.05" --filter "AF <= 0.95" --filter "MAX_AF >= 0.001" --filter "AF >= 0.001"',
        ),
        extra=config.get("germline_vcf", {}).get("extra", ""),
    log:
        "cnv_sv/germline_vcf/{sample}_{type}.germline.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/germline_vcf/{sample}_{type}.germline.vcf.benchmark.tsv",
            config.get("germline_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("germline_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("germline_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("germline_vcf", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("germline_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("germline_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("germline_vcf", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("germline_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/germline_vcf.yaml"
    message:
        "{rule}: Create a germline only vcf cnv_sv/germline_vcf/{wildcards.sample}_{wildcards.type}.germline.vcf"
    shell:
        "(filter_vep -o {output.vcf} {params.filter} {params.extra}) &> {log}"
