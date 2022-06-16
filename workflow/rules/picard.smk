__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule picard_update_vcf_sequence_dictionary:
    input:
        vcf="cnv_sv/pindel/{sample}.noContig.vcf",
        fasta=config["reference"]["fasta"],
    output:
        temp("cnv_sv/pindel/{sample}.vcf"),
    params:
        extra=config.get("picard_update_vcf_sequence_dictionary", {}).get("extra", ""),
    log:
        "cnv_sv/picard/update_vcf_sequence_dictionary/{sample}.output.log",
    benchmark:
        repeat(
            "cnv_sv/picard/update_vcf_sequence_dictionary/{sample}.output.benchmark.tsv",
            config.get("picard_update_vcf_sequence_dictionary", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("picard_update_vcf_sequence_dictionary", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("picard_update_vcf_sequence_dictionary", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("picard_update_vcf_sequence_dictionary", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("picard_update_vcf_sequence_dictionary", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("picard_update_vcf_sequence_dictionary", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("picard_update_vcf_sequence_dictionary", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("picard_update_vcf_sequence_dictionary", {}).get("container", config["default_container"])
    conda:
        "../envs/picard_update_vcf_sequence_dictionary.yaml"
    message:
        "{rule}: Update cnv_sv/pindel/{wildcards.sample}.noContig.vcf to include contigs."
    shell:
        "(picard UpdateVcfSequenceDictionary INPUT={input.vcf} SD={input.fasta} OUTPUT={output}) &> {log}"
