# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule GATK_cnv_callCopyRatioSegments:
    input:
        "cnv/GATK_cnv_modelSegments/{sample}_{type}.clean.cr.seg",
    output:
        temp("cnv/GATK_cnv_callCopyRatioSegments/{sample}_{type}.clean.calledCNVs.seg"),
        temp("cnv/GATK_cnv_callCopyRatioSegments/{sample}_{type}.clean.calledCNVs.igv.seg"),
    params:
        extra=config.get("GATK_cnv_callCopyRatioSegments", {}).get("extra", ""),
    log:
        "cnv/GATK_cnv_callCopyRatioSegments/{sample}_{type}.clean.calledCNVs.seg.log",
    benchmark:
        repeat(
            "cnv/GATK_cnv_callCopyRatioSegments/{sample}_{type}.clean.calledCNVs.seg.benchmark.tsv",
            config.get("GATK_cnv_callCopyRatioSegments", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("GATK_cnv_callCopyRatioSegments", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("GATK_cnv_callCopyRatioSegments", {}).get("container", config["default_container"])
    conda:
        "../envs/GATK_cnv_callCopyRatioSegments.yaml"
    message:
        "{rule}: Use GATK_cnv to obtain cnv/GATK_cnv_callCopyRatioSegments/{wildcards.sample}_{wildcards.type}.clean.calledCNVs.seg"
    shell:
        "(gatk --java-options '-Xmx4g' CallCopyRatioSegments "
        "--input {input} "
        "--output {output} "
        "{params.extra}) &> {log}"
