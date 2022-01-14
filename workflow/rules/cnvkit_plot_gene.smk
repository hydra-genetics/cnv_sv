# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_plot_gene:
    input:
        segments="cnv_sv/cnvkit_call/{sample}/{sample}_{type}.cns",
        segment_regions="cnv_sv/cnvkit_call/{sample}/{sample}_{type}.cnr",
        vcf="cnv_sv/germline_vcf/{sample}_{type}.germline.vcf",
    output:
        segment_regions=temp("cnv_sv/cnvkit_call/{sample}/{sample}_{type}.gene.cnr"),
        plot=temp("cnv_sv/cnvkit_plot_gene/{sample}_{type}_gene.png"),
    params:
        extra=config.get("cnvkit_plot_gene", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_plot_gene/{sample}_{type}.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_plot_gene/{sample}_{type}.benchmark.tsv",
            config.get("cnvkit_plot_gene", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_plot_gene", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("cnvkit_plot_gene", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_plot_gene", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("cnvkit_plot_gene", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_plot_gene", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_plot_gene", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("cnvkit_plot_gene", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit_plot_gene.yaml"
    message:
        "{rule}: Plot cnvs into cnv_sv/cnvkit_plot_gene/{wildcards.sample}_{wildcards.type}.png"
    shell:
        "(sed -e 's/_\(Exon\|Intron\|Intergenic\)[0-9]*_[a-zA-Z0-9_.]*//g' {input.segment_regions} > {output.segment_regions} && "
        "cnvkit.py scatter {output.segment_regions} -s {input.segments} -v {input.vcf} -o {output.plot} -g {input.gene_list} {params.extra}) &> {log}"
