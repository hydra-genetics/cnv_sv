__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2025, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule melt_melt:
    input:
        input1="...",
    output:
        output1="cnv_sv/melt_melt/{sample}_{type}.output.txt",
    params:
        extra=config.get("melt_melt", {}).get("extra", ""),
    log:
        "cnv_sv/melt_melt/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/melt_melt/{sample}_{type}.output.benchmark.tsv",
            config.get("melt_melt", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("melt_melt", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("melt_melt", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("melt_melt", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("melt_melt", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("melt_melt", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("melt_melt", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("melt_melt", {}).get("container", config["default_container"])
    message:
        "{rule}: do stuff on {input.input1}"
    wrapper:
        "..."
