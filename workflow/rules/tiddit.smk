__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2022, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule tiddit:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
    output:
        dir=directory("cnv_sv/tiddit/{sample}_{type}"),
        vcf="cnv_sv/tiddit/{sample}_{type}.vcf"
    params:
        extra=config.get("tiddit", {}).get("extra", ""),
    log:
        "cnv_sv/tiddit/{sample}_{type}.log"
    benchmark:
        repeat("cnv_sv/tiddit/{sample}_{type}.output.benchmark.tsv",
        config.get("tiddit", {}).get("benchmark_repeats", 1))
    threads: config.get("tiddit", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("tiddit", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("tiddit", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("tiddit", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("tiddit", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("tiddit", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("tiddit", {}).get("container", config["default_container"])
    conda:
        "../envs/tiddit.yaml"
    message:
       "{rule}: Run tiddit on {wildcards.sample}_{wildcards.type}"
    shell:
        "TIDDIT.py --sv --bam {input.bam} -o {output.dir} &> {log}"""
