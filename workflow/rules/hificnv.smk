__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2025, Andrei Guliaev"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"


rule hificnv:
    input:
        bam=lambda wildcards: get_input_bam(wildcards, "T")[0],
    output:
        vcf=temp("cnv_sv/hificnv/{sample}_{type}.vcf.gz"),
        bw=temp("cnv_sv/hificnv/{sample}_{type}.depth.bw"),
        bedgraph=temp("cnv_sv/hificnv/{sample}_{type}.copynum.bedgraph"),
    params:
        ref=config.get("reference", {}).get("fasta", ""),
        exclude=config.get("hificnv", {}).get("exclude", ""),
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
        "OUTPUT={output.vcf} && "
        "PREFIX=${{OUTPUT%.vcf.gz}} && "
        "hificnv --bam {input.bam} "
        "--ref {params.ref} "
        "--threads {threads} "
        "--exclude {params.exclude} "
        "--output-prefix $PREFIX "
        "&> {log.call} && "
        "mv $PREFIX.{wildcards.sample}'.vcf.gz' $PREFIX'.vcf.gz' &> {log.mv_vcf} && "
        "mv $PREFIX.{wildcards.sample}'.depth.bw' $PREFIX'.depth.bw' &> {log.mv_bw} && "
        "mv $PREFIX.{wildcards.sample}'.copynum.bedgraph' $PREFIX'.copynum.bedgraph' &> {log.mv_bedgraph}"
