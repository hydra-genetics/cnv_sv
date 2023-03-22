__author__ = "Padraic Corcoran, Jessika Nordin"
__copyright__ = "Copyright 2023, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule upd:
    input:
        father=""
        index=""
        mother=""
        vcf="{trio_id}.vep_annotated.vcf.gz",
    output:
        bed=temp("{trio_id}.upd_regions.bed"),
    params:
        extra=config.get("upd", {}).get("extra", ""),
    log:
        "{file}.upd.log",
    benchmark:
        repeat("{file}.upd.tsv", config.get("upd", {}).get("benchmark_repeats", 1))
    threads: config.get("upd", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("upd", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("upd", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("upd", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("upd", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("upd", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("upd", {}).get("container", config["default_container"])
    conda:
        "../envs/upd.yaml"
    message:
        "{rule}: Use a trio vcf {input.vcf} to find uniparental disomy"
    shell:
        "(upd "
        "--vcf {input.vcf} "
        "--proband {input.index} "
        "--mother {input.mother} "
        "--father {input.father} "
        "{params.extra} "
        "regions "
        "--out {output.bed}) "
        "&> {log}"
