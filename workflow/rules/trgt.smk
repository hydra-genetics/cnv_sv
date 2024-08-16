__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2024, Padraic Corcoran"
__email__ = "padraic.corocran@scilifelab.uu.se"
__license__ = "GPL-3"


rule trgt_genotype:
    input:
        bam="alignment/minimap2/{sample}_{type}.bam",
        bai="alignment/minimap2/{sample}_{type}.bam.bai",
        bed=config.get("trgt_genotype", {}).get("bed", ""),
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf="cnv_sv/trgt_genotype/{sample}_{type}.vcf.gz",
        bam="cnv_sv/trgt_genotype/{sample}_{type}.spanning.bam",
    params:
        extra=config.get("trgt_genotype", {}).get("extra", ""),
        karyotype=get_karyotype,
        prefix="cnv_sv/trgt_genotype/{sample}_{type}",
        sample_name="{sample}_{type}",
    log:
        "cnv_sv/trgt_genotype/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/trgt_genotype/{sample}_{type}.vcf.benchmark.tsv", config.get("trgt_genotype", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("trgt_genotype", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("trgt_genotype", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("trgt_genotype", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("trgt_genotype", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("trgt_genotype", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("trgt_genotype", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("trgt_genotype", {}).get("container", config["default_container"])
    message:
        "{rule}: Run trgt genotype on {input.bam}"
    shell:
        "trgt genotype "
        "--genome {input.ref} "
        "--repeats {input.bed} "
        "--reads {input.bam} "
        "--karyotype {params.karyotype} "
        "--sample-name {params.sample_name} "
        "{params.extra} "
        "--output-prefix {params.prefix} "
        "--threads {threads} &> {log}"


rule trgt_bam_sort:
    input:
        bam="cnv_sv/trgt_genotype/{sample}_{type}.spanning.bam",
    output:
        bam="cnv_sv/trgt_genotype/{sample}_{type}.spanning.sorted.bam",
        idx="cnv_sv/trgt_genotype/{sample}_{type}.spanning.sorted.bam.bai",
    params:
        extra=config.get("trgt_bam_sort", {}).get("extra", ""),
    log:
        "cnv_sv/trgt_genotype/{sample}_{type}.spanning.sorted.bam.log",
    benchmark:
        repeat(
            "cnv_sv/trgt_genotype/{sample}_{type}.output.benchmark.tsv",
            config.get("trgt_bam_sort", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("trgt_bam_sort", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("trgt_bam_sort", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("trgt_bam_sort", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("trgt_bam_sort", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("trgt_bam_sort", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("trgt_bam_sort", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("trgt_bam_sort", {}).get("container", config["default_container"])
    message:
        "{rule}: Sort and index {input.bam} with samtools"
    wrapper:
        "v3.10.2/bio/samtools/sort"


rule trgt_plot:
    input:
        bam="cnv_sv/trgt_genotype/{sample}_{type}.spanning.sorted.bam",
        bai="cnv_sv/trgt_genotype/{sample}_{type}.spanning.sorted.bam.bai",
        bed=config.get("trgt_genotype", {}).get("bed", ""),
        ref=config.get("reference", {}).get("fasta", ""),
        vcf="cnv_sv/trgt_genotype/{sample}_{type}.vcf.gz",
    output:
        image=expand("cnv_sv/trgt_plot/{{sample}}_{{type}}_{{locus}}.{ext}", ext=config.get("trgt_plot", {}).get("image", "svg")),
    params:
        extra=config.get("trgt_plot", {}).get("extra", ""),
        image_type=config.get("trgt_plot", {}).get("image", "svg"),
        plot_type=config.get("trgt_plot", {}).get("plot_type", "allele"),
        show=config.get("trgt_plot", {}).get("show", "motifs"),
    log:
        "cnv_sv/trgt_plot/{sample}_{type}_{locus}.output.log",
    benchmark:
        repeat(
            "cnv_sv/trgt_plot/{sample}_{type}_{locus}.output.benchmark.tsv",
            config.get("trgt_plot", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("trgt_plot", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("trgt_plot", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("trgt_plot", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("trgt_plot", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("trgt_plot", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("trgt_plot", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("trgt_plot", {}).get("container", config["default_container"])
    message:
        "{rule}: Run trgt plot on {input.bam}"
    shell:
        "trgt plot --genome {input.ref} "
        "--repeats {input.bed} "
        "--spanning-reads {input.bam} "
        "--vcf {input.vcf} "
        "--repeat-id {wildcards.locus} "
        "--plot-type {params.plot_type} "
        "--image {output.image} "
        "--show {params.show} "
        "{params.extra} &> {log} || true "
        "&& touch {output.image}"

