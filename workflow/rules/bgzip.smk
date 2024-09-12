__author__ = "Pádraic Corcoran"
__copyright__ = "Copyright 2024, Pádraic Corcoran"
__email__ = "padraic.corcoran@sciliflab.uu.se"
__license__ = "GPL-3"


rule bgzip:
    input:
        vcf="{file}.vcf",
    output:
        gz=temp("{file}.vcf.gz"),
    log:
        "{file}.vcf.gz.log",
    benchmark:
        repeat("{file}.vcf.gz.benchmark.tsv", config.get("bgzip", {}).get("benchmark_repeats", 1))
    threads: config.get("bgzip", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bgzip", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bgzip", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bgzip", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bgzip", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bgzip", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bgzip", {}).get("container", config["default_container"])
    message:
        "{rule}: bgzip {input.vcf}"
    wrapper:
        "v1.3.1/bio/bgzip"
