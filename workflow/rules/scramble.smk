__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2025, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule scramble_scramble:
    input:
        input1="...",
    output:
        output1="cnv_sv/scramble_scramble/{sample}_{type}.output.txt",
    params:
        extra=config.get("scramble_scramble", {}).get("extra", ""),
    log:
        "cnv_sv/scramble_scramble/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/scramble_scramble/{sample}_{type}.output.benchmark.tsv",
            config.get("scramble_scramble", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("scramble_scramble", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("scramble_scramble", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("scramble_scramble", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("scramble_scramble", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("scramble_scramble", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("scramble_scramble", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("scramble_scramble", {}).get("container", config["default_container"])
    message:
        "{rule}: do stuff on {input.input1}"
    wrapper:
        "..."
