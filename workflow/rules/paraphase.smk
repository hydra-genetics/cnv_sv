__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2024, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule paraphase:
    input:
        bam=lambda wildcards: get_input_aligned_bam(wildcards, config)[0],
        bai=lambda wildcards: get_input_aligned_bam(wildcards, config)[1],
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        bam="cnv_sv/paraphase/paraphase_{sample}_{type}/{sample}_{type}.paraphase.bam",
        bai="cnv_sv/paraphase/paraphase_{sample}_{type}/{sample}_{type}.paraphase.bam.bai",
        json="cnv_sv/paraphase/paraphase_{sample}_{type}/{sample}_{type}.paraphase.json",
        vcf=expand(
            "cnv_sv/paraphase/paraphase_{{sample}}_{{type}}/{{sample}}_{{type}}_paraphase_vcfs/{{sample}}_{{type}}_{gene}.vcf",
            gene=config.get("paraphase", {}).get("genes", ""),
        ),
    params:
        extra=config.get("paraphase", {}).get("extra", ""),
        genome=config.get("paraphase", {}).get("genome", "38"),
        prefix=lambda wildcards, output: "{}_{}".format(wildcards.sample, wildcards.type),
        out=lambda wildcards, output: os.path.dirname(output.bam),
        gene=lambda wildcards: "--gene " + ",".join(config.get("paraphase", {}).get("genes", "")),
    log:
        "cnv_sv/paraphase/{sample}_{type}.output.log",
    benchmark:
        repeat("cnv_sv/paraphase/{sample}_{type}.output.benchmark.tsv", config.get("paraphase", {}).get("benchmark_repeats", 1))
    threads: config.get("paraphase", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("paraphase", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("paraphase", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("paraphase", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("paraphase", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("paraphase", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("paraphase", {}).get("container", config["default_container"])
    message:
        "{rule}: run paraphase on  {input.bam}"
    shell:
        "paraphase "
        "--bam {input.bam} "
        "--reference {input.ref} "
        "--prefix {params.prefix} "
        "--out {params.out} "
        "{params.gene} "
        "{params.extra} "
        "--threads {threads} "
        "--genome {params.genome} &> {log}"
