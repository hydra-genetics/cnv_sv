__author__ = "Arielle R Munters, Padraic Corcoran"
__copyright__ = "Copyright 2023, Arielle R Munters, Padraic Corcoran"
__email__ = "arielle.munters@scilifelab.uu.se, padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule exomedepth_call:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bedfile=config.get("exomedepth_call", {}).get("bedfile", ""),
        sex="qc/peddy/peddy.sex_check.csv",
    output:
        exon=temp("cnv_sv/exomedepth_call/{sample}_{type}.RData"),
        txt=temp("cnv_sv/exomedepth_call/{sample}_{type}.txt"),
    params:
        extra=config.get("exomedepth_call", {}).get("extra", ""),
        ref_count=lambda wildcards, input: get_exomedepth_ref(wildcards, input.sex),
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
    message:
        "{rule}: run exomedepth cnv_sv/{rule}/{wildcards.sample}_{wildcards.type}.input"
    script:
        "../scripts/exomedepth.R"
