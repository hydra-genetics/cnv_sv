__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2024, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule juli_call:
    input:
        input1="...",
    output:
        output1="cnv_sv/juli_call/{sample}_{type}.output.txt",
    params:
        extra=config.get("juli_call", {}).get("extra", ""),
    log:
        "cnv_sv/juli_call/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/juli_call/{sample}_{type}.output.benchmark.tsv",
            config.get("juli_call", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("juli_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("juli_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("juli_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("juli_call", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("juli_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("juli_call", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("juli_call", {}).get("container", config["default_container"])
    message:
        "{rule}: do stuff on {input.input1}"
    wrapper:
        "..."


rule juli_annotate:
    input:
        input1="...",
    output:
        output1="cnv_sv/juli_annotate/{sample}_{type}.output.txt",
    params:
        extra=config.get("juli_annotate", {}).get("extra", ""),
    log:
        "cnv_sv/juli_annotate/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/juli_annotate/{sample}_{type}.output.benchmark.tsv",
            config.get("juli_annotate", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("juli_annotate", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("juli_annotate", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("juli_annotate", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("juli_annotate", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("juli_annotate", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("juli_annotate", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("juli_annotate", {}).get("container", config["default_container"])
    message:
        "{rule}: do stuff on {input.input1}"
    wrapper:
        "..."


rule juli_filter:
    input:
        input1="...",
    output:
        output1="cnv_sv/juli_filter/{sample}_{type}.output.txt",
    params:
        extra=config.get("juli_filter", {}).get("extra", ""),
    log:
        "cnv_sv/juli_filter/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/juli_filter/{sample}_{type}.output.benchmark.tsv",
            config.get("juli_filter", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("juli_filter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("juli_filter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("juli_filter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("juli_filter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("juli_filter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("juli_filter", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("juli_filter", {}).get("container", config["default_container"])
    message:
        "{rule}: do stuff on {input.input1}"
    wrapper:
        "..."
