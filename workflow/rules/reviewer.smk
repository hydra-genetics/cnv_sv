__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule reviewer_genrate_locus_list:
    input:
        cat=config.get("expansionhunter", {}).get("variant_catalog", ""),
        vcf="cnv_sv/expansionhunter/{sample}_{type}.vcf",
    output:
        txt=temp("cnv_sv/reviewer/{sample}_{type}_locus_list.txt"),
    log:
        "cnv_sv/reviewer/{sample}_{type}_locus_list.txt.log",
    benchmark:
        repeat(
            "cnv_sv/reviewer/{sample}_{type}_locus_list.txt.benchmark.tsv",
            config.get("reviewer_generate_locus_list", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("reviewer_generate_locus_list", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("reviewer_generate_locus_list", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("reviewer_generate_locus_list", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("reviewer_generate_locus_list", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("reviewer_generate_locus_list", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("reviewer_generate_locus_list", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("reviewer_generate_locus_list", {}).get("container", config["default_container"])
    conda:
        "../envs/reviewer_generate_locus_list.yaml"
    message:
        "{rule}: Generate a locus list for REViewer from {input.vcf}"
    script:
        "../scripts/reviewer_generate_locus_list.py"


rule reviewer:
    input:
        bam="cnv_sv/expansionhunter/{sample}_{type}_realigned.sorted.bam",
        bai="cnv_sv/expansionhunter/{sample}_{type}_realigned.sorted.bam.bai",
        cat=config.get("expansionhunter", {}).get("variant_catalog", ""),
        loci="cnv_sv/reviewer/{sample}_{type}_locus_list.txt",
        ref=config.get("reference", {}).get("fasta", ""),
        vcf="cnv_sv/expansionhunter/{sample}_{type}.vcf",
    output:
        temp(directory("cnv_sv/reviewer/{sample}_{type}/")),
        temp("cnv_sv/reviewer/{sample}_{type}/{sample}_{type}.metrics.tsv"),
        temp("cnv_sv/reviewer/{sample}_{type}/{sample}_{type}.phasing.tsv"),
    params:
        extra=config.get("reviewer", {}).get("extra", ""),
        in_locus=lambda wildcards, input: get_locus_str(input.loci),
        prefix=lambda wildcards, output: "{}/{}_{}/{}_{}".format(
            os.path.split(output[0])[0], wildcards.sample, wildcards.type, wildcards.sample, wildcards.type
        ),
    log:
        "cnv_sv/reviewer/{sample}_{type}/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/reviewer/{sample}_{type}/{sample}_{type}.output.benchmark.tsv",
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
