# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule gatk_cnv_modelSegments:
    input:
        denoisedCopyRatio="cnv/gatk_cnv_denoiseReadCounts/{sample}_{type}.clean.denoisedCR.tsv",
        allelicCounts="cnv/gatk_cnv_collectAllelicCounts/{sample}_{type}.clean.allelicCounts.tsv",
    output:
        temp("cnv/gatk_cnv_modelSegments/{sample}_{type}.clean.modelFinal.seg"),
        temp("cnv/gatk_cnv_modelSegments/{sample}_{type}.clean.cr.seg"),
        temp("cnv/gatk_cnv_modelSegments/{sample}_{type}.clean.af.igv.seg"),
        temp("cnv/gatk_cnv_modelSegments/{sample}_{type}.clean.cr.igv.seg"),
        temp("cnv/gatk_cnv_modelSegments/{sample}_{type}.clean.hets.tsv"),
        temp("cnv/gatk_cnv_modelSegments/{sample}_{type}.clean.modelBegin.cr.param"),
        temp("cnv/gatk_cnv_modelSegments/{sample}_{type}.clean.modelBegin.af.param"),
        temp("cnv/gatk_cnv_modelSegments/{sample}_{type}.clean.modelBegin.seg"),
        temp("cnv/gatk_cnv_modelSegments/{sample}_{type}.modelFinal.af.param"),
        temp("cnv/gatk_cnv_modelSegments/{sample}_{type}.modelFinal.cr.param"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        outprefix="{sample}_{type}.clean",
        extra=config.get("gatk_cnv_modelSegments", {}).get("extra", ""),
    log:
        "cnv/gatk_cnv_modelSegments/{sample}_{type}.clean.modelFinal.seg.log",
    benchmark:
        repeat(
            "cnv/gatk_cnv_modelSegments/{sample}_{type}.clean.modelFinal.seg.benchmark.tsv",
            config.get("gatk_cnv_modelSegments", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_cnv_modelSegments", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("gatk_cnv_modelSegments", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk_cnv_modelSegments.yaml"
    message:
        "{rule}: Use gatk_cnv to obtain cnv/gatk_cnv_modelSegments/{wildcards.sample}_{wildcards.type}.clean.modelFinal.seg"
    shell:
        "(gatk --java-options '-Xmx4g' ModelSegments "
        "--denoised-copy-ratios {input.denoisedCopyRatio} "
        "--allelic-counts {input.allelicCounts} "
        "--output {params.outdir} "
        "--output-prefix {params.outprefix}"
        "{params.extra}) &> {log}"
