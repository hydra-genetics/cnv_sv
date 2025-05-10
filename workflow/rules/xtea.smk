__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2025, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule xtea_xtea:
    input:
        input1="...",
    output:
        output1="cnv_sv/xtea_xtea/{sample}_{type}.output.txt",
    params:
        extra=config.get("xtea_xtea", {}).get("extra", ""),
    log:
        "cnv_sv/xtea_xtea/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/xtea_xtea/{sample}_{type}.output.benchmark.tsv",
            config.get("xtea_xtea", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("xtea_xtea", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("xtea_xtea", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("xtea_xtea", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("xtea_xtea", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("xtea_xtea", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("xtea_xtea", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("xtea_xtea", {}).get("container", config["default_container"])
    message:
        "{rule}: do stuff on {input.input1}"
    wrapper:
        "..."
