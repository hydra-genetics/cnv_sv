# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


if config.get("cnvkit_call", {}).get("method", "hybrid") == "hybrid":

    rule cnvkit_call:
        input:
            bam="alignment/bwa_mem/{sample}_{type}.bam",
            bai="alignment/bwa_mem/{sample}_{type}.bam.bai",
            cnv_reference=config["cnvkit_call"]["normal_reference"],
        output:
            regions=temp("cnv/cnvkit_call/{sample}/{sample}_{type}.cnr"),
            segments=temp("cnv/cnvkit_call/{sample}/{sample}_{type}.cns"),
            segments_called=temp("cnv/cnvkit_call/{sample}/{sample}_{type}.call.cns"),
            bins=temp("cnv/cnvkit_call/{sample}/{sample}_{type}.bintest.cns"),
            target_coverage=temp("cnv/cnvkit_call/{sample}/{sample}_{type}.targetcoverage.cnn"),
            antitarget_coverage=temp("cnv/cnvkit_call/{sample}/{sample}_{type}.antitargetcoverage.cnn"),
        params:
            outdir=lambda wildcards, output: os.path.dirname(output[0]),
            extra=config.get("cnvkit_call", {}).get("extra", ""),
        log:
            "cnv/cnvkit_call/{sample}/{sample}_{type}.log",
        benchmark:
            repeat(
                "cnv/cnvkit_call/{sample}/{sample}_{type}.benchmark.tsv",
                config.get("cnvkit_call", {}).get("benchmark_repeats", 1),
            )
        threads: config.get("cnvkit_call", config["default_resources"]).get("threads", config["default_resources"]["threads"])
        container:
            config.get("cnvkit_call", {}).get("container", config["default_container"])
        conda:
            "../envs/cnvkit_call.yaml"
        message:
            "{rule}: Use cnvkit to call cnvs in {wildcards.sample}/{wildcards.sample}_{wildcards.type}"
        shell:
            "(cnvkit.py batch {input.bam} "
            "-r {input.cnv_reference} "
            "-d {params.outdir} "
            "{params.extra}) &> {log}"
