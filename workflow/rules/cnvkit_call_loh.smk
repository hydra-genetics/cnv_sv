# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_call_loh:
    input:
        segment="cnv/cnvkit_call/{sample}/{sample}_{type}.cns",
        vcf="cnv/germline_vcf/{sample}_{type}.germline.vcf",
    output:
        segment="cnv/cnvkit_call_loh/{sample}_{type}.loh.cns",
    params:
        TC=lambda wildcards: get_sample(samples, wildcards)["TC"],
        extra=config.get("cnvkit_call_loh", {}).get("extra", ""),
    log:
        "cnv/cnvkit_call_loh/{sample}_{type}.loh.cns.log",
    benchmark:
        repeat(
            "cnv/cnvkit_call_loh/{sample}_{type}.loh.cns.benchmark.tsv",
            config.get("cnvkit_call_loh", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_call_loh", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("cnvkit_call_loh", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit_call_loh.yaml"
    message:
        "{rule}: Call cnvs with loh info into cnv/cnvkit_call_loh/{wildcards.sample}_{wildcards.type}.loh.cns"
    shell:
        "(cnvkit.py call {input.segment} -v {input.vcf} -o {output.segment} --purity {params.TC}) &> {log}"
