__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2024, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule sniffles2_call:
    input:
        bam="alignment/minimap2_align/{sample}_{type}.bam",
        bai="alignment/minimap2_align/{sample}_{type}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf="cnv_sv/sniffles2_call/{sample}_{type}.vcf",
        snf="cnv_sv/sniffles2_call/{sample}_{type}.snf"
    params:
        tandem_repeats=get_tr_bed,
        extra=config.get("sniffles2_call", {}).get("extra", ""),
    log:
        "cnv_sv/sniffles2_call/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/sniffles2_call/{sample}_{type}.output.benchmark.tsv",
            config.get("sniffles2_call", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("sniffles2_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sniffles2_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sniffles2_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sniffles2_call", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sniffles2_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sniffles2_call", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sniffles2_call", {}).get("container", config["default_container"])
    message:
        "{rule}: Calls SVs on {input.bam} with sniffles"
    shell:
        "sniffles -i {input.bam} "
        "--reference {input.ref} "
        "-t {threads} "
        "{params.tandem_repeats} "
        "{params.extra} "
        "--vcf {output.vcf} "
        "--snf {output.snf} &> {log}"
