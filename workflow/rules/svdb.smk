# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule svdb_merge:
    input:
        vcfs=get_vcfs_for_svdb_merge,
    output:
        vcf=temp("cnv_sv/svdb_merge/{sample}_{type}.{tc_method}.merged.vcf"),
    params:
        extra=config.get("svdb_merge", {}).get("extra", ""),
        overlap=config.get("svdb_merge", {}).get("overlap", 0.6),
        bnd_distance=config.get("svdb_merge", {}).get("bnd_distance", 10000),
        vcfs=lambda wildards: get_vcfs_for_svdb_merge(wildards, add_suffix=True)
    log:
        "cnv_sv/svdb_merge/{sample}_{type}.{tc_method}.merged.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/svdb_merge/{sample}_{type}.{tc_method}.merged.benchmark.tsv",
            config.get("svdb_merge", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("svdb_merge", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("svdb_merge", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("svdb_merge", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("svdb_merge", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("svdb_merge", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("svdb_merge", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("svdb_merge", {}).get("container", config["default_container"])
    message:
        "{rule}: merges vcf files from different cnv callers into {output.vcf}"
    shell:
        "(svdb --merge "
        "--vcf {params.vcfs} "
        "--bnd_distance {params.bnd_distance} "
        "--overlap {params.overlap} "
        "{params.extra} "
        "> {output.vcf}) 2> {log}"


rule svdb_query:
    input:
        vcf="cnv_sv/svdb_merge/{sample}_{type}.{tc_method}.merged.vcf",
    output:
        vcf=temp("cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.vcf"),
    params:
        db_string=config.get("svdb_query", {}).get("db_string", ""),
        extra=config.get("svdb_query", {}).get("extra", ""),
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0][:-6],
    log:
        "cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.log",
    benchmark:
        repeat(
            "cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.benchmark.tsv",
            config.get("svdb_query", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("svdb_query", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("svdb_query", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("svdb_query", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("svdb_query", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("svdb_query", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("svdb_query", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("svdb_query", {}).get("container", config["default_container"])
    message:
        "{rule}: use svdb database to filter cnvs into {output.vcf}"
    shell:
        "(svdb --query "
        "--query_vcf {input.vcf} "
        "{params.db_string} "
        "--prefix {params.prefix} "
        "{params.extra}) &> {log}"
