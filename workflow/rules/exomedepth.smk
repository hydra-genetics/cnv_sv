__author__ = "Arielle R Munters, Padraic Corcoran"
__copyright__ = "Copyright 2023, Arielle R Munters, Padraic Corcoran"
__email__ = "arielle.munters@scilifelab.uu.se, padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


checkpoint exomedepth_sex:
    input:
        "qc/peddy/peddy.sex_check.csv",
    output:
        "qc/peddy/peddy.sex_check-checkpoint.csv",
    params:
        extra=config.get("exomedepth_sex", {}).get("extra", ""),
    log:
        "qc/peddy/peddy.sex_check-checkpoint.csv.log",
    benchmark:
        repeat(
            "qc/peddy/peddy.sex_check-checkpoint.csv.benchmark.tsv",
            config.get("exomedepth_sex", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("exomedepth_sex", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("exomedepth_sex", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("exomedepth_sex", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("exomedepth_sex", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("exomedepth_sex", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("exomedepth_sex", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("exomedepth_sex", {}).get("container", config["default_container"])
    message:
        "{rule}: checkpoint with {input}"
    shell:
        "cp {input} {output}"


rule exomedepth_call:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bedfile=config.get("exomedepth_call", {}).get("bedfile", ""),
        sex="qc/peddy/peddy.sex_check-checkpoint.csv",
        ref_count=get_exomedepth_ref,
    output:
        exon=temp("cnv_sv/exomedepth_call/{sample}_{type}.RData"),
        txt=temp("cnv_sv/exomedepth_call/{sample}_{type}.txt"),
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
    message:
        "{rule}: run exomedepth_call with {input.bam}, and ref count {input.ref_count}"
    script:
        "../scripts/exomedepth.R"
