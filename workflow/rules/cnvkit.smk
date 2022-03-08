# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf, Martin Rippin"
__copyright__ = "Copyright 2022, Jonas Almlöf, Martin Rippin"
__email__ = "jonas.almlof@scilifelab.uu.se, martin.rippin@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_batch:
    input:
        bam="alignment/merge_bam/{sample}_{type}.bam",
        bai="alignment/merge_bam/{sample}_{type}.bam.bai",
        cnv_reference=config.get("cnvkit_batch", {}).get("normal_reference", ""),
    output:
        regions=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr"),
        segments=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns"),
        segments_called=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.call.cns"),
        bins=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.bintest.cns"),
        target_coverage=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.targetcoverage.cnn"),
        antitarget_coverage=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.antitargetcoverage.cnn"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        method=config.get("cnvkit_batch", {}).get("method", "hybrid"),
        extra=config.get("cnvkit_batch", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.benchmark.tsv",
            config.get("cnvkit_batch", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_batch", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("cnvkit_batch", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_batch", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("cnvkit_batch", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_batch", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_batch", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("cnvkit_batch", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit.yaml"
    message:
        "{rule}: Use cnvkit to call cnvs in {wildcards.sample}/{wildcards.sample}_{wildcards.type}"
    shell:
        "(cnvkit.py batch {input.bam} "
        "-r {input.cnv_reference} "
        "-d {params.outdir} "
        "-m {params.method} "
        "{params.extra}) &> {log}"


rule cnvkit_call:
    input:
        segment="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
        vcf="cnv_sv/germline_vcf/{sample}_{type}.germline.vcf",
    output:
        segment=temp("cnv_sv/cnvkit_call/{sample}_{type}.loh.cns"),
    params:
        TC=lambda wildcards: get_sample(samples, wildcards)["tumor_content"],
        extra=config.get("cnvkit_call", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_call/{sample}_{type}.loh.cns.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_call/{sample}_{type}.loh.cns.benchmark.tsv",
            config.get("cnvkit_call", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("cnvkit_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_call", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("cnvkit_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_call", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("cnvkit_call", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit.yaml"
    message:
        "{rule}: Call cnvs with loh info into cnv_sv/cnvkit_call/{wildcards.sample}_{wildcards.type}.loh.cns"
    shell:
        "(cnvkit.py call {input.segment} "
        "-v {input.vcf} "
        "-o {output.segment} "
        "--purity {params.TC} "
        "{params.extra}) &> {log}"


rule cnvkit_diagram:
    input:
        cns="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
        cnr="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr",
    output:
        pdf=temp("cnv_sv/cnvkit_diagram/{sample}_{type}.pdf"),
    params:
        extra=config.get("cnvkit_diagram", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_diagram/{sample}_{type}.pdf.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_diagram/{sample}_{type}.pdf.benchmark.tsv",
            config.get("cnvkit_diagram", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_diagram", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("cnvkit_diagram", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_diagram", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("cnvkit_diagram", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_diagram", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_diagram", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("cnvkit_diagram", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit.yaml"
    message:
        "{rule}: chromosome plot cnv_sv/cnvkit_scatter/{wildcards.sample}_{wildcards.type}.pdf"
    shell:
        "(cnvkit.py diagram {input.cnr} "
        "-s {input.cns} "
        "-o {output.pdf} "
        "{params.extra}) &> {log}"


rule cnvkit_scatter:
    input:
        segments="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
        segment_regions="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr",
        vcf="cnv_sv/germline_vcf/{sample}_{type}.germline.vcf",
    output:
        plot=temp("cnv_sv/cnvkit_scatter/{sample}_{type}.png"),
    params:
        extra=config.get("cnvkit_scatter", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_scatter/{sample}_{type}.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_scatter/{sample}_{type}.benchmark.tsv",
            config.get("cnvkit_scatter", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_scatter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("cnvkit_scatter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_scatter", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("cnvkit_scatter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_scatter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_scatter", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("cnvkit_scatter", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit.yaml"
    message:
        "{rule}: Plot cnvs into cnv_sv/cnvkit_scatter/{wildcards.sample}_{wildcards.type}.png"
    shell:
        "(cnvkit.py scatter {input.segment_regions} "
        "-s {input.segments} "
        "-v {input.vcf} "
        "-o {output.plot} "
        "{params.extra}) &> {log}"
