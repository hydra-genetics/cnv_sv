
__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"

rule expansionhunter:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
        cat=config.get("expansionhunter", {}).get("variant_catalog", ""),
    output:
        vcf="cnv_sv/expansionhunter/{sample}_{type}.vcf",
        json="cnv_sv/expansionhunter/{sample}_{type}.json",
        bam=temp("cnv_sv/expansionhunter/{sample}_{type}_realigned.bam"),
    params:
        prefix = lambda wildcards, output: os.path.split(output.vcf)[0],
        extra=config.get("expansionhunter", {}).get("extra", ""),
    log:
        "cnv_sv/expansionhunter/{sample}_{type}.output.log"
    benchmark:
        repeat("cnv_sv/expansionhunter/{sample}_{type}.output.benchmark.tsv",
        config.get("expansionhunter", {}).get("benchmark_repeats", 1))
    threads: config.get("expansionhunter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("expansionhunter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("expansionhunter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("expansionhunter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("expansionhunter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("expansionhunter", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("expansionhunter", {}).get("container", config["default_container"])
    conda:
        "../envs/expansionhunter.yaml"
    message:
       "{rule}: Run ExpansionHunter on {wildcards.sample}_{wildcards.type}"
    shell:
        "ExpansionHunter --reads {input.bam} "
        "--reference {input.ref} "
        "--threads {threads} "
        "--variant-catalog {input.cat} "
        "--output-prefix {params.prefix}/{wildcards.sample}_{wildcards.type} &> {log}"


rule samtools_sort:
    input:
        "cnv_sv/expansionhunter/{sample}_{type}_realigned.bam",
    output:
        temp("cnv_sv/expansionhunter/{sample}_{type}_realigned.sorted.bam"),
    params:
        extra=config.get("samtools_sort", {}).get("extra", ""),
    log:
        "{path_file}.bam.sort.log",
    benchmark:
        repeat(
            "{path_file}.bam.sort.benchmark.tsv",
            config.get("samtools_sort", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("samtools_sort", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("samtools_sort", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("samtools_sort", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("samtools_sort", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("samtools_sort", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("samtools_sort", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("samtools_sort", {}).get("container", config["default_container"])
    conda:
        "../envs/samtools.yaml"
    message:
        "{rule}: sort bam file {input} using samtools"
    wrapper:
        "v1.3.2/bio/samtools/sort"

rule samtools_index:
    input:
        "cnv_sv/expansionhunter/{sample}_{type}_realigned.sorted.bam",
    output:
        "cnv_sv/expansionhunter/{sample}_{type}_realigned.sorted.bam.bai",
    params:
        extra=config.get("samtools_index", {}).get("extra", ""),
    log:
        "{file}.bam.bai.log",
    benchmark:
        repeat(
            "{file}.bam.bai.benchmark.tsv",
            config.get("samtools_index", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("samtools_index", {}).get("container", config["default_container"])
    threads: config.get("samtools_index", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("samtools_index", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("samtools_index", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("samtools_index", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("samtools_index", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("samtools_index", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../envs/samtools.yaml"
    message:
        "{rule}: create index for cnv_sv/expansionhunter/{sample}_{type}_realigned.sorted.bam"
    wrapper:
        "v1.1.0/bio/samtools/index"

rule reviewer:
    input:
        bam="cnv_sv/expansionhunter/{sample}_{type}_realigned.sorted.bam",
        vcf="cnv_sv/expansionhunter/{sample}_{type}.vcf"
        ref=config.get("reference", {}).get("fasta", ""),
        cat=config.get("reviewer", {}).get("variant_catalog", ""),
    output:
        vcf="cnv_sv/reviewer/{sample}_{type}.vcf",
        json="cnv_sv/reviewer/{sample}_{type}.json",
        bam="cnv_sv/reviewer/{sample}_{type}_realigned.bam",
    params:
        prefix = lambda wildcards, output: os.path.split(output.vcf)[0],
        extra=config.get("reviewer", {}).get("extra", ""),
        in_locus=get_locus_list
    log:
        "cnv_sv/reviewer/{sample}_{type}.output.log"
    benchmark:
        repeat("cnv_sv/reviewer/{sample}_{type}.output.benchmark.tsv",
        config.get("reviewer", {}).get("benchmark_repeats", 1))
    threads: config.get("reviewer", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("reviewer", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("reviewer", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("reviewer", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("reviewer", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("reviewer", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("reviewer", {}).get("container", config["default_container"])
    conda:
        "../envs/reviewer.yaml"
    message:
       "{rule}: Run reviewer on {wildcards.sample}_{wildcards.type}"
    shell:
        "reviewer --reads {input.bam} "
        "--reference {input.ref} "
        "--catalog {input.cat} "
        "--locus {parmas.in_locus}"
        "--output-prefix {params.prefix}/{wildcards.sample}_{wildcards.type} &> {log}"
