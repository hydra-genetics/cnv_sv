# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule GATK_cnv_denoiseReadCounts:
    input:
        hdf5PoN=config["GATK_cnv_denoiseReadCounts"]["normal_reference"],
        hdf5Tumor="cnv/GATK_cnv_collectReadCounts/{sample}_{type}.counts.hdf5",
    output:
        denoisedCopyRatio=temp("cnv/GATK_cnv_denoiseReadCounts/{sample}_{type}.clean.denoisedCR.tsv"),
        stdCopyRatio=temp("cnv/GATK_cnv_denoiseReadCounts/{sample}_{type}.clean.standardizedCR.tsv"),
    params:
        extra=config.get("GATK_cnv_denoiseReadCounts", {}).get("extra", ""),
    log:
        "cnv/GATK_cnv_denoiseReadCounts/{sample}_{type}.clean.denoisedCR.tsv.log",
    benchmark:
        repeat(
            "cnv/GATK_cnv_denoiseReadCounts/{sample}_{type}.clean.denoisedCR.tsv.benchmark.tsv",
            config.get("GATK_cnv_denoiseReadCounts", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("GATK_cnv_denoiseReadCounts", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("GATK_cnv_denoiseReadCounts", {}).get("container", config["default_container"])
    conda:
        "../envs/GATK_cnv_denoiseReadCounts.yaml"
    message:
        "{rule}: Use GATK_cnv to obtain cnv/GATK_cnv_denoiseReadCounts/{wildcards.sample}_{wildcards.type}.clean.denoisedCR.tsv"
    shell:
        "(gatk --java-options '-Xmx4g' DenoiseReadCounts -I {input.hdf5Tumor} "
        "--count-panel-of-normals {input.hdf5PoN} "
        "--standardized-copy-ratios {output.stdCopyRatio} "
        "--denoised-copy-ratios {output.denoisedCopyRatio} "
        "{params.exra}) &> {log}"
