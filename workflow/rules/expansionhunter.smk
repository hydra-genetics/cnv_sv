
__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"

rule expansionhunter:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
        cat=config.get("expansionhunter", {}).get("variant_catalog", ""),
    output:
        vcf="cnv_sv/expansionhunter/{sample}_{type}.vcf",
        json="cnv_sv/expansionhunter/{sample}_{type}.json",
        bam="cnv_sv/expansionhunter/{sample}_{type}_realigned.bam",
    params:
        prefix = lambda wildcards, output: os.path.split(output.vcf)[0],
        extra=config.get("expansionhunter", {}).get("extra", ""),
    log:
        "cnv_sv/expansionhunter/{sample}_{type}.output.log"
    benchmark:
        repeat("cnv_sv/expansionhunter/{sample}_{type}.output.benchmark.tsv",
        config.get("expansionhunter", {}).get("benchmark_repeats", 1))
    threads: config.get("expansionhunter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("expansionhunter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("expansionhunter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("expansionhunter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("expansionhunter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("expansionhunter", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("expansionhunter", {}).get("container", config["default_container"])
    conda:
        "../envs/expansionhunter.yaml"
    message:
       "{rule}: Run ExpansionHunter on {wildcards.sample}_{wildcards.type}"
    shell:
        "ExpansionHunter --reads {input.bam} "
        "--reference {input.ref} "
        "--threads {threads} "
        "--variant-catalog {input.cat} "
        "--output-prefix {params.prefix}/{wildcards.sample}_{wildcards.type} &> {log}"
