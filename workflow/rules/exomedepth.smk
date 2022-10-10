__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule exomedepth_call:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bedfile=config.get("exomedepth_call", {}).get("bedfile", ""),
        ref_count=config.get("exomedepth_call", {}).get("ref_count", ""),
    output:
        aed=temp("cnv_sv/exomedepth_call/{sample}_{type}.aed"),
        aggregated_result=temp("cnv_sv/exomedepth_call/{sample}_{type}.SV.txt"),
        exon=temp("cnv_sv/exomedepth_call/{sample}_{type}.RData"),
        result=temp("cnv_sv/exomedepth_call/{sample}_{type}.txt"),
    params:
        extra=config.get("exomedepth_call", {}).get("extra", ""),
    log:
        "cnv_sv/exomedepth_call/{sample}_{type}.txt.log",
    benchmark:
        repeat(
            "cnv_sv/exomedepth_call/{sample}_{type}.txt.benchmark.tsv",
            config.get("exomedepth_call", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("exomedepth_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("exomedepth_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("exomedepth_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("exomedepth_call", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("exomedepth_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("exomedepth_call", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("exomedepth_call", {}).get("container", config["default_container"])
    conda:
        "../envs/exomedepth.yaml"
    message:
        "{rule}: run exomedepth cnv_sv/{rule}/{wildcards.sample}_{wildcards.type}.input"
    script:
        "../scripts/exomedepth.R"


rule exomedepth_ed_filter:
    input:
        conifer=config.get("exomedepth_ed_filter", {}).get("conifer", ""),
        ED_common=config.get("exomedepth_ed_filter", {}).get("ED_common", ""),
        exon="cnv_sv/exomedepth_call/{sample}_{type}.RData",
        WW25=config.get("exomedepth_ed_filter", {}).get("WW25", ""),
    output:
        aed=temp("cnv_sv/exomedepth_ed_filter/{sample}_{type}.filter.aed"),
        aggregated_result=temp("cnv_sv/exomedepth_ed_filter/{sample}_{type}.SV.filter.txt"),
    params:
        extra=config.get("exomedepth_ed_filter", {}).get("extra", ""),
    log:
        "cnv_sv/exomedepth_ed_filter/{sample}_{type}.SV.filter.txt.log",
    benchmark:
        repeat(
            "cnv_sv/exomedepth_ed_filter/{sample}_{type}.SV.filter.txt.benchmark.tsv",
            config.get("exomedepth_ed_filter", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("exomedepth_ed_filter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("exomedepth_ed_filter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("exomedepth_ed_filter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("exomedepth_ed_filter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("exomedepth_ed_filter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("exomedepth_ed_filter", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("exomedepth_ed_filter", {}).get("container", config["default_container"])
    conda:
        "../envs/exomedepth.yaml"
    message:
        "{rule}: filter exomedepth cnv_sv/{rule}/{wildcards.sample}_{wildcards.type}.RData"
    script:
        "../scripts/exomedepth_filter.R"
