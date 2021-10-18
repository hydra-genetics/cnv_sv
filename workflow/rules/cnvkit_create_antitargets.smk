# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_create_antitargets:
    input:
        bed="cnv/cnvkit_create_targets/cnvkit.target.bed",
    output:
        bed=temp("cnv/cnvkit_create_antitargets/cnvkit.antitarget.bed"),
    log:
        "cnv/cnvkit_create_antitargets/cnvkit.antitarget.bed.log",
    benchmark:
        repeat(
            "cnv/cnvkit_create_antitargets/cnvkit.antitarget.bed.benchmark.tsv",
            config.get("cnvkit_create_antitargets", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_create_antitargets", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("cnvkit_create_antitargets", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit_create_antitargets.yaml"
    message:
        "{rule}: Use cnvkit to create antitarget bedfile cnvkit.antitarget.bed"
    shell:
        "(cnvkit.py antitarget {input.bed} -o {output.bed}) &> {log}"
