__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2025, Andrei Guliaev"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"


rule hificnv:
    input:
        bam=lambda wildcards: get_input_bam(wildcards)[0],
        bai=lambda wildcards: get_input_bam(wildcards)[1],
    output:
        vcf=temp("cnv_sv/hificnv/{sample}_{type}.vcf.gz"),
        bw=temp("cnv_sv/hificnv/{sample}_{type}.depth.bw"),
        bedgraph=temp("cnv_sv/hificnv/{sample}_{type}.copynum.bedgraph"),
    params:
        ref=config.get("reference", {}).get("fasta", ""),
        exclude=config.get("hificnv", {}).get("exclude", ""),
        output_prefix= lambda wildcards,output: str(output.vcf).replace(".vcf.gz",""),
    log:
        call="cnv_sv/hificnv/{sample}_{type}.vcf.gz.log",
        mv_vcf="cnv_sv/hificnv/{sample}_{type}.vcf.gz.mv.log",
        mv_bw="cnv_sv/hificnv/{sample}_{type}.depth.bw.mv.log",
        mv_bedgraph="cnv_sv/hificnv/{sample}_{type}.copynum.bedgraph.mv.log",
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
        "{rule}: Calculating copy number variants on {input.bam} with HiFiCNV"
    shell:
        "hificnv --bam {input.bam} "
        "--ref {params.ref} "
        "--threads {threads} "
        "--exclude {params.exclude} "
        "--output-prefix {params.output_prefix} "
        "&> {log.call} && "
        "mv {params.output_prefix}.{wildcards.sample}_{wildcards.type}.vcf.gz {params.output_prefix}.vcf.gz &> {log.mv_vcf} && "
        "mv {params.output_prefix}.{wildcards.sample}_{wildcards.type}.depth.bw {params.output_prefix}.depth.bw &> {log.mv_bw} && "
        "mv {params.output_prefix}.{wildcards.sample}_{wildcards.type}.copynum.bedgraph {params.output_prefix}.copynum.bedgraph &> {log.mv_bedgraph}"
