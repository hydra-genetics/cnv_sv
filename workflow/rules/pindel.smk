__author__ = "Martin Rippin"
__copyright__ = "Copyright 2022, Martin Rippin"
__email__ = "martin.rippin@scilifelab.uu.se"
__license__ = "GPL-3"


rule pindel_call:
    input:
        config="cnv_sv/pindel/{sample}.cfg",
        bam="alignment/merge_bam/{sample}_T.bam",
        bai="alignment/merge_bam/{sample}_T.bam.bai",
        ref=config["reference"]["fasta"],
    output:
        pindel=temp(
            expand(
                "cnv_sv/pindel/{{sample}}_{ext}",
                ext=[
                    "BP",
                    "CloseEndMapped",
                    "D",
                    "INT_final",
                    "INV",
                    "LI",
                    "RP",
                    "SI",
                    "TD",
                ],
            )
        ),
    params:
        prefix=lambda wildcards: "cnv_sv/pindel/%s" % (wildcards.sample),
        extra=config.get("pindel", {}).get("extra", ""),
    log:
        "cnv_sv/pindel/{sample}_pindel_call.log",
    benchmark:
        repeat(
            "cnv_sv/pindel/{sample}_pindel_call.benchmark.tsv",
            config.get("pindel_call", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("pindel_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("pindel_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pindel_call", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("pindel_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pindel_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pindel_call", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("pindel_call", {}).get("container", config["default_container"])
    conda:
        "../envs/pindel.yaml"
    message:
        "{rule}: Detect breakpoints in {wildcards.sample}"
    wrapper:
        "v1.2.0/bio/pindel/call"


rule pindel2vcf:
    input:
        pindel=expand(
            "cnv_sv/pindel/{{sample}}_{ext}",
            ext=[
                "BP",
                "CloseEndMapped",
                "D",
                "INT_final",
                "INV",
                "LI",
                "RP",
                "SI",
                "TD",
            ],
        ),
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("cnv_sv/pindel/{sample}.vcf"),
    params:
        refname=config.get("pindel2vcf", {}).get("refname", "'Genome Reference Consortium Human Build 38'"),
        refdate=config.get("pindel2vcf", {}).get("refdate", "20131217"),
        extra=config.get("pindel2vcf", {}).get("extra", ""),
    log:
        "cnv_sv/pindel/{sample}_pindel2vcf.log",
    benchmark:
        repeat(
            "cnv_sv/pindel/{sample}_pindel2vcf.benchmark.tsv",
            config.get("pindel2vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("pindel2vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("pindel2vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pindel2vcf", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("pindel2vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pindel2vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pindel2vcf", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("pindel2vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/pindel.yaml"
    message:
        "{rule}: Convert pindel output to vcf for {wildcards.sample}"
    wrapper:
        "v1.2.0/bio/pindel/pindel2vcf"
