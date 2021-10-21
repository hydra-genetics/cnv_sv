# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule GATK_cnv_modelSegments:
    input:
        denoisedCopyRatio="cnv/GATK_cnv_modelSegments/{sample}_{type}.clean.denoisedCR.tsv",
        allelicCounts="cnv/GATK_cnv_modelSegments/{sample}_{type}.clean.allelicCounts.tsv",
    output:
        temp("cnv/GATK_cnv_modelSegments/{sample}_{type}.clean.modelFinal.seg"),
        temp("cnv/GATK_cnv_modelSegments/{sample}_{type}.clean.sr.seg"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        outPrefix="{sample}_{type}.clean",
        extra=config.get("GATK_cnv_modelSegments", {}).get("extra", ""),
    log:
        "cnv/GATK_cnv_modelSegments/{sample}_{type}.clean.modelFinal.seg.log",
    benchmark:
        repeat(
            "cnv/GATK_cnv_modelSegments/{sample}_{type}.clean.modelFinal.seg.benchmark.tsv",
            config.get("GATK_cnv_modelSegments", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("GATK_cnv_modelSegments", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("GATK_cnv_modelSegments", {}).get("container", config["default_container"])
    conda:
        "../envs/GATK_cnv_modelSegments.yaml"
    message:
        "{rule}: Use GATK_cnv to obtain cnv/GATK_cnv_modelSegments/{wildcards.sample}_{wildcards.type}.clean.modelFinal.seg"
    shell:
        "(gatk --java-options '-Xmx4g' ModelSegments "
        "--denoised-copy-ratios {input.denoisedCopyRatio} "
        "--allelic-counts {input.allelicCounts} "
        "--output {params.outDir} "
        "--output-prefix {params.outPrefix}"
        "{params.extra}) &> {log}"
