__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2022, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvpytor_readdepth:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
    output:
        pytor="cnv_sv/cnvpytor/{sample}_{type}.pytor",
    params:
        extra=config.get("cnvpytor", {}).get("extra", ""),
    log:
        "cnv_sv/cnvpytor/{sample}_{type}.log"
    benchmark:
        repeat("cnv_sv/cnvpytor/{sample}_{type}.output.benchmark.tsv",
        config.get("cnvpytor", {}).get("benchmark_repeats", 1))
    threads: config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvpytor", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvpytor", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvpytor", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvpytor", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvpytor", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvpytor.yaml"
    message:
       "{rule}: Run cnvpytor on {wildcards.sample}_{wildcards.type}"
    shell:
        "cnvpytor -root {output.pytor} -rd {input.bam} &> {log}"
