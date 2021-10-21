# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_coverage:
    input:
        bed="cnv/cnvkit_create_{target}argets/{target}arget.bed",
        bam="alignment/bwa_mem/{sample}_{type}.bam",
        bai="alignment/bwa_mem/{sample}_{type}.bam.bai",
    output:
        cnn=temp("cnv/cnvkit_coverage/{sample}_{type}.{target}argetcoverage.cnn"),
    log:
        "cnv/cnvkit_coverage/{sample}_{type}.{target}argetcoverage.cnn.log",
    benchmark:
        repeat(
            "cnv/cnvkit_coverage/{sample}_{type}.{target}argetcoverage.cnn.benchmark.tsv",
            config.get("cnvkit_coverage", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_coverage", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("cnvkit_coverage", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit_coverage.yaml"
    message:
        "{rule}: Use cnvkit to calculate coverage in {wildcards.sample}_{wildcards.type}.{wildcards.target}argetcoverage.cnn"
    shell:
        "(cnvkit.py coverage {input.bam} {input.bed} -o {output.cnn}) &> {log}"
