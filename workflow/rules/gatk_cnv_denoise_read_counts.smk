# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule gatk_cnv_denoise_read_counts:
    input:
        hdf5PoN=config["gatk_cnv_denoise_read_counts"]["normal_reference"],
        hdf5Tumor="cnv_sv/gatk_cnv_collect_read_counts/{sample}_{type}.counts.hdf5",
    output:
        denoisedCopyRatio=temp("cnv_sv/gatk_cnv_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv"),
        stdCopyRatio=temp("cnv_sv/gatk_cnv_denoise_read_counts/{sample}_{type}.clean.standardizedCR.tsv"),
    params:
        extra=config.get("gatk_cnv_denoise_read_counts", {}).get("extra", ""),
    log:
        "cnv_sv/gatk_cnv_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv.log",
    benchmark:
        repeat(
            "cnv_sv/gatk_cnv_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv.benchmark.tsv",
            config.get("gatk_cnv_denoise_read_counts", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_cnv_denoise_read_counts", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("gatk_cnv_denoise_read_counts", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_cnv_denoise_read_counts", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("gatk_cnv_denoise_read_counts", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_cnv_denoise_read_counts", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_cnv_denoise_read_counts", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("gatk_cnv_denoise_read_counts", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk_cnv_denoise_read_counts.yaml"
    message:
        "{rule}: Use gatk_cnv to obtain cnv_sv/gatk_cnv_denoise_read_counts/{wildcards.sample}_{wildcards.type}.clean.denoisedCR.tsv"
    shell:
        "(gatk --java-options '-Xmx4g' DenoiseReadCounts -I {input.hdf5Tumor} "
        "--count-panel-of-normals {input.hdf5PoN} "
        "--standardized-copy-ratios {output.stdCopyRatio} "
        "--denoised-copy-ratios {output.denoisedCopyRatio} "
        "{params.extra}) &> {log}"
