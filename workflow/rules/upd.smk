__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2023, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule bcftools_biallelic:
    input:
        vcf="{file}.vep_annotated.vcf",
    output:
        vcf=temp("{file}.biallelic.vcf.gz"),
    params:
        extra=config.get("bcftools_biallelic", {}).get("extra", ""),
    log:
        "{file}.biallelic.log",
    benchmark:
        repeat("{file}.biallelic.benchmark.tsv", config.get("cnvpytor", {}).get("benchmark_repeats", 1))
    threads: config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvpytor", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvpytor", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvpytor", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvpytor", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvpytor_readdepth", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvpytor.yaml"
    message:
        "{rule}: Get biallelic variants from {input.vcf}"
    shell:
        "(bcftools view --min-alleles 2 --max-alleles 2  "
        "--types snps {input.vcf} --apply-filters PASS "
        "--output-type z --output-file {output.vcf}) &> {log}"


rule bcftools_biallelic:
    input:
        vcf="{file}.vep_annotated.vcf",
    output:
        vcf=temp("{file}.biallelic.vcf.gz"),
    params:
        extra=config.get("bcftools_biallelic", {}).get("extra", ""),
    log:
        "{file}.biallelic.log",
    benchmark:
        repeat("{file}.biallelic.benchmark.tsv", config.get("cnvpytor", {}).get("benchmark_repeats", 1))
    threads: config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvpytor", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvpytor", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvpytor", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvpytor", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvpytor_readdepth", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvpytor.yaml"
    message:
        "{rule}: Get biallelic variants from {input.vcf}"
    shell:
        "(bcftools view --min-alleles 2 --max-alleles 2  "
        "--types snps {input.vcf} --apply-filters PASS "
        "--output-type z --output-file {output.vcf}) &> {log}"
