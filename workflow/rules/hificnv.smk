__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2025, Andrei Guliaev"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"

rule hificnv:
    input:
        bam=lambda wildcards: get_longread_bam(wildcards)[0],
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        prefix="cnv_sv/hificnv/{sample}_{type}",
    params:
        exclude=onfig.get("hificnv", {}).get("exclude", ""),
    log:
        "cnv_sv/hificnv/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/hificnv/{sample}_{type}.output.benchmark.tsv",
            config.get("hificnv", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("hificnv", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("hificnv", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("hificnv", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("hificnv", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("hificnv", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("hificnv", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("hificnv", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculates CNVs on {input.bam} with HiFiCNV"
    shell:
        "hificnv --bam {input.bam} "
        "--reference {input.ref} "
        "--threads {threads} "
        "--exclude {params.exclude} "
        "--output_prefix {output.prefix} "
        "&> {log}"