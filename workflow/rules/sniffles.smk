__copyright__ = "Copyright 2023, Nina Hollfelder"
__email__ = "nina_hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


rule sniffles:
    input:
        bam="alignment/minimap2/{sample}_{type}.bam",
        bai="alignment/minimap2/{sample}_{type}.bam.bai",
        fasta=config.get("reference", {}).get("fasta", ""),
    output:
        vcf="cnv_sv/sniffles/{sample}_{type}.vcf.gz",
        snf="cnv_sv/sniffles/{sample}_{type}.snf",
    params:
        extra=config.get("sniffles", {}).get("extra", ""),
    log:
        "cnv_sv/sniffles/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/sniffles/{sample}_{type}.vcf.benchmark.tsv",
            config.get("sniffles", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("sniffles", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sniffles", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sniffles", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sniffles", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sniffles", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sniffles", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sniffles", {}).get("container", config["default_container"])
    message:
        "{rule}: Calls SVs on {input.bam} with sniffles"
    shell:
        "sniffles -i {input.bam} "
        "--reference {input.fasta} "
        "-t {threads} "
        "{params.extra} "
        "--vcf {output.vcf} "
        "--snf {output.snf} &> {log}"
