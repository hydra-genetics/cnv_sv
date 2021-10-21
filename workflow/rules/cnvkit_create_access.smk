# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_create_access:
    input:
        fasta=config["reference"]["fasta"],
    output:
        bed=temp("cnv/cnvkit_create_access/access_excludes.bed"),
    params:
        extra=config.get("cnvkit_create_access", {}).get("extra", ""),
    log:
        "cnv/cnvkit_create_access/access_excludes.bed.log",
    benchmark:
        repeat(
            "cnv/cnvkit_create_access/access_excludes.bed.benchmark.tsv",
            config.get("cnvkit_create_access", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_create_access", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("cnvkit_create_access", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit_create_access.yaml"
    message:
        "{rule}: Use cnvkit to create a bedfile with regions that are inaccessible in sequencing in access_excludes.bed"
    shell:
        "(cnvkit.py access {input.fasta} -o {output.bed}) &> {log}"
