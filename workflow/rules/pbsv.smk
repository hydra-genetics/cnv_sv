__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2025, Andrei Guliaev"
__email__ = "andrei.guliae@scilifelab.uu.se"
__license__ = "GPL-3"


rule pbsv_discover:
    input:
        bam=lambda wildcards: get_input_bam(wildcards, "T")[0],
        bai=lambda wildcards: get_input_bam(wildcards, "T")[1],
    output:
        svsig=temp("cnv_sv/pbsv_discover/{sample}_{type}.svsig.gz"),
    params:
        trf=config.get("pbsv_discover", {}).get("trf", ""),
        extra=config.get("pbsv_discover", {}).get("extra", ""),
    log:
        discover="cnv_sv/pbsv_discover/{sample}_{type}.svsig.gz.log",
    benchmark:
        repeat(
            "cnv_sv/pbsv_discover/{sample}_{type}.output.benchmark.tsv",
            config.get("pbsv_discover", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("pbsv_discover", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("pbsv_discover", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pbsv_discover", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pbsv_discover", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("pbsv_discover", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pbsv_discover", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("pbsv_discover", {}).get("container", config["default_container"])
    message:
        "{rule}: discover SV signatures with pbsv on {input.bam}"
    shell:
        "pbsv discover "
        "--tandem-repeats {params.trf} "
        "--log-file {log.discover} "
        "{input.bam} {output.svsig}"


rule pbsv_call:
    input:
        svsig="cnv_sv/pbsv_discover/{sample}_{type}.svsig.gz",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf=temp("cnv_sv/pbsv_call/{sample}_{type}.vcf"),
    params:
        extra=config.get("pbsv_call", {}).get("extra", ""),
    log:
        call="cnv_sv/pbsv_call/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/pbsv_call/{sample}_{type}.output.benchmark.tsv",
            config.get("pbsv_call", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("pbsv_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("pbsv_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pbsv_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pbsv_call", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("pbsv_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pbsv_call", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("pbsv_call", {}).get("container", config["default_container"])
    message:
        "{rule}: call SV with pbsv from signatures file {input.svsig}"
    shell:
        "pbsv call -j {threads} "
        "--log-file {log.call} "
        "--hifi {input.ref} "
        "{input.svsig} "
        "{output.vcf}"
