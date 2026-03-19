__author__ = "Julia Höglund"
__copyright__ = "Copyright 2025, Julia Höglund"
__email__ = "julia.hoglund@scilifelab.uu.se"
__license__ = "GPL-3"


rule scramble_cluster_identifier:
    input:
        bam=lambda wildcards: get_input_aligned_bam(wildcards, config)[0],
        bai=lambda wildcards: get_input_aligned_bam(wildcards, config)[1],
    output:
        clusters="cnv_sv/scramble_cluster_identifier/{sample}_{type}.clusters.txt",
    params:
        extra=config.get("scramble_cluster_identifier", {}).get("extra", ""),
    log:
        "cnv_sv/scramble_cluster_identifier/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/scramble_cluster_identifier/{sample}_{type}.output.benchmark.tsv",
            config.get("scramble_cluster_identifier", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("scramble_cluster_identifier", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("scramble_cluster_identifier", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("scramble_cluster_identifier", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("scramble_cluster_identifier", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("scramble_cluster_identifier", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("scramble_cluster_identifier", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("scramble_cluster_identifier", {}).get("container", config["default_container"])
    message:
        "{rule}: identify read clusters in {input.bam} with SCRAMble"
    shell:
        "cluster_identifier "
        "{params.extra} "
        "{input.bam} "
        "> {output.clusters} "
        "2> {log}"


rule scramble_cluster_analysis:
    input:
        clusters="cnv_sv/scramble_cluster_identifier/{sample}_{type}.clusters.txt",
    output:
        meis="cnv_sv/scramble_cluster_analysis/{sample}_{type}_MEIs.txt",
    params:
        extra=config.get("scramble_cluster_analysis", {}).get("extra", ""),
        install_dir=config.get("scramble_cluster_analysis", {}).get("install_dir", ""),
        mei_refs=config.get("scramble_cluster_analysis", {}).get("mei_refs", ""),
        out_name=lambda wildcards: f"cnv_sv/scramble_cluster_analysis/{wildcards.sample}_{wildcards.type}",
    log:
        "cnv_sv/scramble_cluster_analysis/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/scramble_cluster_analysis/{sample}_{type}.output.benchmark.tsv",
            config.get("scramble_cluster_analysis", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("scramble_cluster_analysis", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("scramble_cluster_analysis", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("scramble_cluster_analysis", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("scramble_cluster_analysis", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("scramble_cluster_analysis", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("scramble_cluster_analysis", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("scramble_cluster_analysis", {}).get("container", config["default_container"])
    message:
        "{rule}: analyze read clusters from {input.clusters} using SCRAMble"
    shell:
        "WORKDIR=$(pwd); "
        "Rscript --vanilla /usr/local/bin/SCRAMble.R "
        "--out-name $WORKDIR/{params.out_name} "
        "--cluster-file $WORKDIR/{input.clusters} "
        "--install-dir {params.install_dir} "
        "--mei-refs {params.mei_refs} "
        "--eval-meis "
        "{params.extra} "
        "> {log} 2>&1"


rule scramble_vcf:
    input:
        meis="cnv_sv/scramble_cluster_analysis/{sample}_{type}_MEIs.txt",
    output:
        vcf=temp("cnv_sv/scramble_vcf/{sample}_{type}.vcf"),
    params:
        sample_name=lambda wildcards: f"{wildcards.sample}_{wildcards.type}",
        min_alignment_score=config.get("scramble_vcf", {}).get("min_alignment_score", 70),
        min_clipped_reads=config.get("scramble_vcf", {}).get("min_clipped_reads", 5),
        min_alignment_percent_length=config.get("scramble_vcf", {}).get("min_alignment_percent_length", 90),
        min_alignment_percent_identity=config.get("scramble_vcf", {}).get("min_alignment_percent_identity", 90),
        alu_size=config.get("scramble_vcf", {}).get("alu_size", 282),
        sva_size=config.get("scramble_vcf", {}).get("sva_size", 1362),
        l1_size=config.get("scramble_vcf", {}).get("l1_size", 6023),
        cluster_distance=config.get("scramble_vcf", {}).get("cluster_distance", 300),
    log:
        "cnv_sv/scramble_vcf/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/scramble_vcf/{sample}_{type}.vcf.benchmark.tsv",
            config.get("scramble_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("scramble_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("scramble_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("scramble_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("scramble_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("scramble_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("scramble_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("scramble_vcf", {}).get("container", config["default_container"])
    message:
        "{rule}: convert {input.meis} to VCF format"
    script:
        "../scripts/scramble_vcf.py"


rule scramble_sort:
    input:
        vcf="cnv_sv/scramble_vcf/{sample}_{type}.vcf.gz",
        tbi="cnv_sv/scramble_vcf/{sample}_{type}.vcf.gz.tbi",
        fai=config.get("reference", {}).get("fai", ""),
    output:
        vcf=temp("cnv_sv/scramble_vcf/{sample}_{type}.sorted.vcf.gz"),
    log:
        "cnv_sv/scramble_vcf/{sample}_{type}.sorted.vcf.gz.log",
    benchmark:
        repeat(
            "cnv_sv/scramble_vcf/{sample}_{type}.sorted.vcf.gz.benchmark.tsv",
            config.get("scramble_sort", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("scramble_sort", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("scramble_sort", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("scramble_sort", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("scramble_sort", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("scramble_sort", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("scramble_sort", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("scramble_sort", {}).get("container", config["default_container"])
    message:
        "{rule}: sort and reheader {input.vcf} with bcftools"
    shell:
        """
        bcftools reheader \\
        -f {input.fai}  {input.vcf} \\
        | bcftools sort -Oz -o {output.vcf} \\
        &> {log}
        """
