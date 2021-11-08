# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule gatk_cnv_model_segments:
    input:
        denoisedCopyRatio="cnv/gatk_cnv_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv",
        allelicCounts="cnv/gatk_cnv_collect_allelic_counts/{sample}_{type}.clean.allelicCounts.tsv",
    output:
        temp("cnv/gatk_cnv_model_segments/{sample}_{type}.clean.modelFinal.seg"),
        temp("cnv/gatk_cnv_model_segments/{sample}_{type}.clean.cr.seg"),
        temp("cnv/gatk_cnv_model_segments/{sample}_{type}.clean.af.igv.seg"),
        temp("cnv/gatk_cnv_model_segments/{sample}_{type}.clean.cr.igv.seg"),
        temp("cnv/gatk_cnv_model_segments/{sample}_{type}.clean.hets.tsv"),
        temp("cnv/gatk_cnv_model_segments/{sample}_{type}.clean.modelBegin.cr.param"),
        temp("cnv/gatk_cnv_model_segments/{sample}_{type}.clean.modelBegin.af.param"),
        temp("cnv/gatk_cnv_model_segments/{sample}_{type}.clean.modelBegin.seg"),
        temp("cnv/gatk_cnv_model_segments/{sample}_{type}.clean.modelFinal.af.param"),
        temp("cnv/gatk_cnv_model_segments/{sample}_{type}.clean.modelFinal.cr.param"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        outprefix="{sample}_{type}.clean",
        extra=config.get("gatk_cnv_model_segments", {}).get("extra", ""),
    log:
        "cnv/gatk_cnv_model_segments/{sample}_{type}.clean.modelFinal.seg.log",
    benchmark:
        repeat(
            "cnv/gatk_cnv_model_segments/{sample}_{type}.clean.modelFinal.seg.benchmark.tsv",
            config.get("gatk_cnv_model_segments", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_cnv_model_segments", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("gatk_cnv_model_segments", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_cnv_model_segments", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_cnv_model_segments", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk_cnv_model_segments.yaml"
    message:
        "{rule}: Use gatk_cnv to obtain cnv/gatk_cnv_model_segments/{wildcards.sample}_{wildcards.type}.clean.modelFinal.seg"
    shell:
        "(gatk --java-options '-Xmx4g' ModelSegments "
        "--denoised-copy-ratios {input.denoisedCopyRatio} "
        "--allelic-counts {input.allelicCounts} "
        "--output {params.outdir} "
        "--output-prefix {params.outprefix}"
        "{params.extra}) &> {log}"
