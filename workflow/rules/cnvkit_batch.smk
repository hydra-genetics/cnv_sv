# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_batch:
    input:
        bam="alignment/merge_bam/{sample}_{type}.bam",
        bai="alignment/merge_bam/{sample}_{type}.bam.bai",
        cnv_reference=config.get("cnvkit_batch", {}).get("normal_reference", ""),
    output:
        regions=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr"),
        segments=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns"),
        segments_called=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.call.cns"),
        bins=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.bintest.cns"),
        target_coverage=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.targetcoverage.cnn"),
        antitarget_coverage=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.antitargetcoverage.cnn"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        method=config.get("cnvkit_batch", {}).get("method", "hybrid"),
        extra=config.get("cnvkit_batch", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.benchmark.tsv",
            config.get("cnvkit_batch", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_batch", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("cnvkit_batch", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_batch", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("cnvkit_batch", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_batch", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_batch", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("cnvkit_batch", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit_batch.yaml"
    message:
        "{rule}: Use cnvkit to call cnvs in {wildcards.sample}/{wildcards.sample}_{wildcards.type}"
    shell:
        "(cnvkit.py batch {input.bam} "
        "-r {input.cnv_reference} "
        "-d {params.outdir} "
        "-m {params.method} "
        "{params.extra}) &> {log}"
