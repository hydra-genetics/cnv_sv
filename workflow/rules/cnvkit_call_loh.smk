# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_call_loh:
    input:
        segment="cnv_sv/cnvkit_call/{sample}/{sample}_{type}.cns",
        vcf="cnv_sv/germline_vcf/{sample}_{type}.germline.vcf",
    output:
        segment=temp("cnv_sv/cnvkit_call_loh/{sample}_{type}.loh.cns"),
    params:
        TC=lambda wildcards: get_sample(samples, wildcards)["TC"],
        extra=config.get("cnvkit_call_loh", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_call_loh/{sample}_{type}.loh.cns.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_call_loh/{sample}_{type}.loh.cns.benchmark.tsv",
            config.get("cnvkit_call_loh", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_call_loh", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("cnvkit_call_loh", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_call_loh", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("cnvkit_call_loh", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_call_loh", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_call_loh", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("cnvkit_call_loh", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit_call_loh.yaml"
    message:
        "{rule}: Call cnvs with loh info into cnv_sv/cnvkit_call_loh/{wildcards.sample}_{wildcards.type}.loh.cns"
    shell:
        "(cnvkit.py call {input.segment} -v {input.vcf} -o {output.segment} --purity {params.TC} {params.extra}) &> {log}"
