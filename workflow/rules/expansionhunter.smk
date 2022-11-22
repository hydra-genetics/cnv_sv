__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule expansionhunter:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        cat=config.get("expansionhunter", {}).get("variant_catalog", ""),
        ref=config.get("reference", {}).get("fasta", ""),
        sex="qc/peddy/peddy.sex_check.csv",
    output:
        bam=temp("cnv_sv/expansionhunter/{sample}_{type}_realigned.bam"),
        json=temp("cnv_sv/expansionhunter/{sample}_{type}.json"),
        vcf=temp("cnv_sv/expansionhunter/{sample}_{type}.vcf"),
    params:
        extra=config.get("expansionhunter", {}).get("extra", ""),
        prefix=lambda wildcards, output: "{}/{}_{}".format(os.path.split(output.vcf)[0], wildcards.sample, wildcards.type),
        sex=lambda wildards, input: get_peddy_sex(wildards, input.sex),
    log:
        "cnv_sv/expansionhunter/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/expansionhunter/{sample}_{type}.output.benchmark.tsv",
            config.get("expansionhunter", {}).get("benchmark_repeats", 1),
        )
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
        "--sex {params.sex} "
        "--variant-catalog {input.cat} "
        "{params.extra} "
        "--output-prefix {params.prefix} &> {log}"


rule generate_reviewer_locus_list:
    input:
        cat=config.get("expansionhunter", {}).get("variant_catalog", ""),
        vcf="cnv_sv/expansionhunter/{sample}_{type}.vcf",
    output:
        txt=temp("cnv_sv/expansionhunter/{sample}_{type}_locus_list.txt"),
    log:
        "cnv_sv/expansionhunter/{sample}_{type}_locus_list.txt.log",
    benchmark:
        repeat(
            "cnv_sv/expansionhunter/{sample}_{type}_locus_list.txt.benchmark.tsv",
            config.get("generate_reviewer_locus_list", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("generate_reviewer_locus_list", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("generate_reviewer_locus_list", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("generate_reviewer_locus_list", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("generate_reviewer_locus_list", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("generate_reviewer_locus_list", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("generate_reviewer_locus_list", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("generate_reviewer_locus_list", {}).get("container", config["default_container"])
    conda:
        "../envs/generate_reviewer_locus_list.yaml"
    message:
        "{rule}: Generate a locus list for REViewer from {input.vcf}"
    script:
        "../scripts/generate_reviewer_locus_list.py"


rule reviewer:
    input:
        bam="cnv_sv/expansionhunter/{sample}_{type}_realigned.sorted.bam",
        bai="cnv_sv/expansionhunter/{sample}_{type}_realigned.sorted.bam.bai",
        cat=config.get("expansionhunter", {}).get("variant_catalog", ""),
        loci="cnv_sv/expansionhunter/{sample}_{type}_locus_list.txt",
        ref=config.get("reference", {}).get("fasta", ""),
        vcf="cnv_sv/expansionhunter/{sample}_{type}.vcf",
    output:
        temp(directory("cnv_sv/expansionhunter/reviewer/{sample}_{type}/")),
        temp("cnv_sv/expansionhunter/reviewer/{sample}_{type}/{sample}_{type}.metrics.tsv"),
        temp("cnv_sv/expansionhunter/reviewer/{sample}_{type}/{sample}_{type}.phasing.tsv"),
    params:
        extra=config.get("reviewer", {}).get("extra", ""),
        in_locus=lambda wildcards, input: get_locus_str(input.loci),
        prefix=lambda wildcards, input: "{}/{}/{}_{}/{}_{}".format(
            os.path.split(input.vcf)[0], "reviewer", wildcards.sample, wildcards.type, wildcards.sample, wildcards.type
        ),
    log:
        "cnv_sv/expansionhunter/reviewer/{sample}_{type}/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/expansionhunter/reviewer/{sample}_{type}/{sample}_{type}.output.benchmark.tsv",
            config.get("reviewer", {}).get("benchmark_repeats", 1),
        )
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
        "REViewer --reads {input.bam} "
        "--vcf {input.vcf} "
        "--reference {input.ref} "
        "--catalog  {input.cat} "
        "--locus {params.in_locus} "
        "--output-prefix {params.prefix} &> {log}"
