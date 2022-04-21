__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


## regel som automatisk genererar bedfile for exomedepth!

rule exomedepth:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bedfile=config["exomedepth"]["bedfile"],
        ref_count=config["exomedepth"]["ref_count"],
    output:
        temp("cnv_sv/exomedepth/{sample}_{type}.txt"),
        temp("cnv_sv/exomedepth/{sample}_{type}.SV.txt"),
        temp("cnv_sv/exomedepth/{sample}_{type}.aed"),
    params:
        extra=config.get("exomedepth", {}).get("extra", ""),
    log:
        "cnv_sv/exomedepth/{sample}_{type}.output.log"
    benchmark:
        repeat(
            "cnv_sv/exomedepth/{sample}_{type}.output.benchmark.tsv", config.get("exomedepth", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("exomedepth", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("exomedepth", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("exomedepth", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("exomedepth", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("exomedepth", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("exomedepth", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("exomedepth", {}).get("container", config["default_container"])
    conda:
        "../envs/exomedepth.yaml"
    message:
       "{rule}: Run exomedepth cnv_sv/{rule}/{wildcards.sample}_{wildcards.type}.input"
    script:
        "../scripts/exomedepth.R"
