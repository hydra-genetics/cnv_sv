__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2025, Andrei Guliaev"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"


rule hificnv:
    input:
        bam=lambda wildcards: get_input_bam(wildcards, "T")[0],
    output:
        vcf="cnv_sv/hificnv/{sample}.vcf.gz",
        bw="cnv_sv/hificnv/{sample}.depth.bw",
        bedgraph="cnv_sv/hificnv/{sample}.copynum.bedgraph",
    params:
        ref=config.get("reference", {}).get("fasta", ""),
        exclude=config.get("hificnv", {}).get("exclude", ""),
    log:
        "cnv_sv/hificnv/{sample}.vcf.gz.log",
    benchmark:
        repeat(
            "cnv_sv/hificnv/{sample}.output.benchmark.tsv",
            config.get("hificnv", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("hificnv", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("hificnv", {}).get(
            "mem_mb", config["default_resources"]["mem_mb"]
        ),
        mem_per_cpu=config.get("hificnv", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("hificnv", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("hificnv", {}).get(
            "threads", config["default_resources"]["threads"]
        ),
        time=config.get("hificnv", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("hificnv", {}).get("container", config["default_container"])
    message:
        "{rule}: Calculating copy number variants on {input.bam} with HiFiCNV"
    shell:
        "OUTNAME={output.vcf} && "
        "hificnv --bam {input.bam} "
        "--ref {params.ref} "
        "--threads {threads} "
        "--exclude {params.exclude} "
        "--output-prefix ${{OUTNAME%.vcf.gz}} "
        "&> {log}"
