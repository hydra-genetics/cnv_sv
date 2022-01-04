__author__ = "Martin Rippin"
__copyright__ = "Copyright 2022, Martin Rippin"
__email__ = "martin.rippin@scilifelab.uu.se"
__license__ = "GPL-3"


rule generate_pindel_config:
    input:
        bam="alignment/merge_bam/{sample}_T.bam",
        metrics="qc/picard_collect_multiple_metrics/{sample}_T.insert_size_metrics",
    output:
        config=temp("cnv_sv/pindel/{sample}.cfg"),
    log:
        "cnv_sv/pindel/{sample}_generate_pindel_config.log",
    benchmark:
        repeat(
            "cnv_sv/pindel/{sample}_generate_pindel_config.benchmark.tsv",
            config.get("generate_pindel_config", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("generate_pindel_config", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("generate_pindel_config", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("generate_pindel_config", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("generate_pindel_config", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("generate_pindel_config", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("generate_pindel_config", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("generate_pindel_config", {}).get("container", config["default_container"])
    conda:
        "../envs/generate_pindel_config.yaml"
    message:
        "{rule}: Produce config for {wildcards.sample}"
    script:
        "../scripts/generate_pindel_config.py"
