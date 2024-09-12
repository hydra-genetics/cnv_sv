__author__ = "Pádraic Corcoran"
__copyright__ = "Copyright 2024, Pádraic Corcoran"
__email__ = "padraic.corcoran@sciliflab.uu.se"
__license__ = "GPL-3"


rule tabix:
    input:
        gz="{file}.vcf.gz",
    output:
        tbi=temp("{file}.vcf.gz.tbi"),
    params:
        extra=config.get("tabix", {}).get("extra", ""),
    log:
        "{file}.vcf.gz.tbi.log",
    benchmark:
        repeat(
            "{file}.vcf.gz.tbi.benchmark.tsv",
            config.get("tabix", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("tabix", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("tabix", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("tabix", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("tabix", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("tabix", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("tabix", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("tabix", {}).get("container", config["default_container"])
    message:
        "{rule}: index {input.gz}"
    wrapper:
        "v1.3.1/bio/tabix"
