__author__ = "Martin Rippin"
__copyright__ = "Copyright 2022, Martin Rippin"
__email__ = "martin.rippin@scilifelab.uu.se"
__license__ = "GPL-3"


rule pindel:
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
        "cnv_sv/pindel/{sample}_pindel.log",
    benchmark:
        repeat(
            "cnv_sv/pindel/{sample}_pindel.benchmark.tsv",
            config.get("pindel", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("pindel", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("pindel", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pindel", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("pindel", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pindel", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pindel", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("pindel", {}).get("container", config["default_container"])
    conda:
        "../envs/pindel.yaml"
    message:
        "{rule}: Detect breakpoints in {wildcards.sample}"
    wrapper:
        "v0.85.1/bio/pindel/call"
