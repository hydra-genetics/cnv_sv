__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2025, Andrei Guliaev"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"


rule severus_t_only:
    input:
        unpack(lambda wildcards:  get_input_haplotagged_bam(wildcards, config, set_type="T")),
        vntr=config.get("severus_t_only", {}).get("vntr", ""),
        pon=config.get("severus_t_only", {}).get("pon", ""),
    output:
        dir=temp(directory("cnv_sv/severus_t_only/{sample}_{type}/out_dir")),
        b_double=temp("cnv_sv/severus_t_only/{sample}_{type}_breakpoint_double.csv"),
        qual=temp("cnv_sv/severus_t_only/{sample}_{type}_read_qual.txt"),
        all_sv_vcf=temp("cnv_sv/severus_t_only/{sample}_{type}_all_sv.vcf"),
        all_b_clusters=temp("cnv_sv/severus_t_only/{sample}_{type}_all_breakpoint_clusters.tsv"),
        all_b_clusters_list=temp("cnv_sv/severus_t_only/{sample}_{type}_all_breakpoint_clusters_list.tsv"),
        somatic_sv_vcf=temp("cnv_sv/severus_t_only/{sample}_{type}_somatic_sv.vcf"),
        somatic_b_clusters=temp("cnv_sv/severus_t_only/{sample}_{type}_somatic_breakpoint_clusters.tsv"),
        somatic_b_clusters_list=temp("cnv_sv/severus_t_only/{sample}_{type}_somatic_breakpoint_clusters_list.tsv"),
    params:
        extra=config.get("severus_t_only", {}).get("extra", ""),
    log:
        "cnv_sv/severus_t_only/{sample}_{type}.severus_t_only.log",
    benchmark:
        repeat(
            "cnv_sv/severus_t_only/{sample}_{type}.severus_t_only.benchmark.tsv",
            config.get("severus_t_only", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("severus_t_only", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("severus_t_only", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("severus_t_only", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("severus_t_only", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("severus_t_only", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("severus_t_only", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("severus_t_only", {}).get("container", config["default_container"])
    message:
        "{rule}: use Severus in tumor only mode to call SV in {wildcards.sample}_{wildcards.type}"
    shell:
        "( severus --target-bam {input.bam_t} "
        "--out-dir {output.dir} "
        "-t {threads} "
        "--vntr-bed {input.vntr} "
        "--PON {input.pon} "
        "{params.extra} && "
        "cp {output.dir}/breakpoints_double.csv {output.b_double} && "
        "cp {output.dir}/read_qual.txt {output.qual} && "
        "cp {output.dir}/all_SVs/severus_all.vcf {output.all_sv_vcf} && "
        "cp {output.dir}/all_SVs/breakpoint_clusters.tsv {output.all_b_clusters} && "
        "cp {output.dir}/all_SVs/breakpoint_clusters_list.tsv {output.all_b_clusters_list} && "
        "cp {output.dir}/somatic_SVs/severus_somatic.vcf {output.somatic_sv_vcf} && "
        "cp {output.dir}/somatic_SVs/breakpoint_clusters.tsv {output.somatic_b_clusters} && "
        "cp {output.dir}/somatic_SVs/breakpoint_clusters_list.tsv {output.somatic_b_clusters_list} ) "
        "&> {log} "


rule severus_tn:
    input:
        unpack(get_tumor_haplotagged_bam),
        unpack(get_normal_haplotagged_bam),
        vntr=config.get("severus_tn", {}).get("vntr", ""),
    output:
        dir=temp(directory("cnv_sv/severus_tn/{sample}_{type}/out_dir")),
        b_double=temp("cnv_sv/severus_tn/{sample}_{type}_breakpoint_double.csv"),
        qual=temp("cnv_sv/severus_tn/{sample}_{type}_read_qual.txt"),
        all_sv_vcf=temp("cnv_sv/severus_tn/{sample}_{type}_all_sv.vcf"),
        all_b_clusters=temp("cnv_sv/severus_tn/{sample}_{type}_all_breakpoint_clusters.tsv"),
        all_b_clusters_list=temp("cnv_sv/severus_tn/{sample}_{type}_all_breakpoint_clusters_list.tsv"),
        somatic_sv_vcf=temp("cnv_sv/severus_tn/{sample}_{type}_somatic_sv.vcf"),
        somatic_b_clusters=temp("cnv_sv/severus_tn/{sample}_{type}_somatic_breakpoint_clusters.tsv"),
        somatic_b_clusters_list=temp("cnv_sv/severus_tn/{sample}_{type}_somatic_breakpoint_clusters_list.tsv"),
    params:
        extra=config.get("severus_tn", {}).get("extra", ""),
    log:
        "cnv_sv/severus_tn/{sample}_{type}.severus_tn.log",
    benchmark:
        repeat(
            "cnv_sv/severus_tn/{sample}_{type}.severus_tn.benchmark.tsv",
            config.get("severus_tn", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("severus_tn", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("severus_tn", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("severus_tn", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("severus_tn", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("severus_tn", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("severus_tn", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("severus_tn", {}).get("container", config["default_container"])
    message:
        "{rule}: use Severus in tumor-normal mode to call SV in {wildcards.sample}_{wildcards.type}"
    shell:
        "( severus --target-bam {input.bam_t} "
        "--control-bam {input.bam_n} "
        "--out-dir {output.dir} "
        "-t {threads} "
        "--vntr-bed {input.vntr} "
        "{params.extra} && "
        "cp {output.dir}/breakpoints_double.csv {output.b_double} && "
        "cp {output.dir}/read_qual.txt {output.qual} && "
        "cp {output.dir}/all_SVs/severus_all.vcf {output.all_sv_vcf} && "
        "cp {output.dir}/all_SVs/breakpoint_clusters.tsv {output.all_b_clusters} && "
        "cp {output.dir}/all_SVs/breakpoint_clusters_list.tsv {output.all_b_clusters_list} && "
        "cp {output.dir}/somatic_SVs/severus_somatic.vcf {output.somatic_sv_vcf} && "
        "cp {output.dir}/somatic_SVs/breakpoint_clusters.tsv {output.somatic_b_clusters} && "
        "cp {output.dir}/somatic_SVs/breakpoint_clusters_list.tsv {output.somatic_b_clusters_list} ) "
        "&> {log} "
