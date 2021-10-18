# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_create_targets:
    input:
        bed=config["reference"]["design_bedfile"],
    output:
        bed=temp("cnv/cnvkit_create_targets/target.bed"),
    log:
        "cnv/cnvkit_create_targets/target.bed.log",
    benchmark:
        repeat(
            "cnv/cnvkit_create_targets/target.bed.benchmark.tsv",
            config.get("cnvkit_create_targets", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_create_targets", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("cnvkit_create_targets", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit_create_targets.yaml"
    message:
        "{rule}: Use cnvkit to create target bedfile target.bed"
    shell:
        "(cnvkit.py target --split {input.bed} -o {output.bed}) &> {log}"
