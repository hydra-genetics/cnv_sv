__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2025, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule scanitd:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf="cnv_sv/scanitd/{sample}_{type}.vcf",
    params:
        region_bed=config.get("scanitd", {}).get("region_bed", ""),
        extra=config.get("scanitd", {}).get("extra", ""),
    log:
        "cnv_sv/scanitd/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/scanitd/{sample}_{type}.vcf.benchmark.tsv",
            config.get("scanitd", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("scanitd", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("scanitd", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("scanitd", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("scanitd", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("scanitd", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("scanitd", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("scanitd", {}).get("container", config["default_container"])
    message:
        "{rule}: call ITD in supplied regions into {output.vcf}"
    shell:
        "(scanitd "
        "-i {input.bam} "
        "-r {input.ref} "
        "-o {output.vcf} "
        "-t {params.region_bed}) "
        "&> {log}"