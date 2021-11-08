# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule germline_vcf:
    input:
        vcf="snv_indel/ensemble_vcf/{sample}_{type}.ensembled.vcf.gz",
    output:
        vcf=temp("cnv/germline_vcf/{sample}_{type}.germline.vcf"),
    params:
        filter=config.get("germline_vcf", {}).get(
            "filter",
            '--filter "DP > 50" --filter "AF >= 0.05" --filter "AF <= 0.95" --filter "MAX_AF >= 0.001" --filter "AF >= 0.001"',
        ),
        extra=config.get("germline_vcf", {}).get("extra", ""),
    log:
        "cnv/germline_vcf/{sample}_{type}.germline.vcf.log",
    benchmark:
        repeat(
            "cnv/germline_vcf/{sample}_{type}.germline.vcf.benchmark.tsv",
            config.get("germline_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("germline_vcf", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("germline_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/germline_vcf.yaml"
    message:
        "{rule}: Create a germline only vcf cnv/germline_vcf/{wildcards.sample}_{wildcards.type}.germline.vcf"
    shell:
        "(zcat {input.vcf} | "
        "filter_vep -o {output.vcf} {params.filter} {params.extra}) &> {log}"
