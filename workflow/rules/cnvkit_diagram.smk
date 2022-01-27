__author__ = "Martin Rippin"
__copyright__ = "Copyright 2022, Martin Rippin"
__email__ = "martin.rippin@scilifelab.uu.se"
__license__ = "GPL-3"


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
        "../envs/cnvkit_diagram.yaml"
    message:
        "{rule}: chromosome plot cnv_sv/cnvkit_scatter/{wildcards.sample}_{wildcards.type}.pdf"
    shell:
        "(cnvkit.py diagram {input.cnr} -s {input.cns} -o {output.pdf} {params.extra}) &> {log}"
