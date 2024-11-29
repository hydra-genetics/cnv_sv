__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2024, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule sawfish_discover:
    input:
        bam=lambda wildcards: get_longread_bam(wildcards)[0],
        bai=lambda wildcards: get_longread_bam(wildcards)[1],
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        temp("cnv_sv/sawfish_discover/{sample}_{type}/assembly.regions.bed"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/candidate.sv.bcf"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/candidate.sv.bcf.csi"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/contig.alignment.bam"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/contig.alignment.bam.csi"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/debug.breakpoint_clusters.bed"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/debug.cluster.refinement.txt"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/depth.bw"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/depth.mpack"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/discover.settings.json"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/genome.gclevels.mpack"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/max.depth.bed"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/run.stats.json"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/sample.gcbias.mpack"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/sawfish.log"),
        temp("cnv_sv/sawfish_discover/{sample}_{type}/expected.copy.number.bed")
        if config.get("sawfish_discover", {}).get("expected_cn", False)
        else [],
    params:
        extra=config.get("sawfish_discover", {}).get("extra", ""),
        expected_cn=get_expected_cn,
        out_dir="cnv_sv/sawfish_discover/{sample}_{type}",
    log:
        "cnv_sv/sawfish_discover/{sample}_{type}.sawfish_discover.log",
    benchmark:
        repeat(
            "cnv_sv/sawfish_discover/{sample}_{type}.sawfish_discover.benchmark.tsv",
            config.get("sawfish_discover", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("sawfish_discover", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sawfish_discover", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sawfish_discover", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sawfish_discover", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sawfish_discover", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sawfish_discover", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sawfish_discover", {}).get("container", config["default_container"])
    message:
        "{rule}: Run sawfish discover {input.bam}"
    shell:
        "sawfish discover "
        "{params.extra} "
        "--threads {threads} "
        "--ref {input.ref} "
        "--bam {input.bam} "
        "--output-dir {params.out_dir} "
        "{params.expected_cn} &> {log}"


rule sawfish_joint_call:
    input:
        bam="alignment/pbmm2_align/{sample}_{type}.bam",
        bai="alignment/pbmm2_align/{sample}_{type}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
        asm_bed="cnv_sv/sawfish_discover/{sample}_{type}/assembly.regions.bed",
        bcf="cnv_sv/sawfish_discover/{sample}_{type}/candidate.sv.bcf",
        csi="cnv_sv/sawfish_discover/{sample}_{type}/candidate.sv.bcf.csi",
        contig_bam="cnv_sv/sawfish_discover/{sample}_{type}/contig.alignment.bam",
        contig_csi="cnv_sv/sawfish_discover/{sample}_{type}/contig.alignment.bam.csi",
        cluster_bed="cnv_sv/sawfish_discover/{sample}_{type}/debug.breakpoint_clusters.bed",
        refine_txt="cnv_sv/sawfish_discover/{sample}_{type}/debug.cluster.refinement.txt",
        dp_bw="cnv_sv/sawfish_discover/{sample}_{type}/depth.bw",
        dp_mpack="cnv_sv/sawfish_discover/{sample}_{type}/depth.mpack",
        settings_json="cnv_sv/sawfish_discover/{sample}_{type}/discover.settings.json",
        cn_bed=(
            "cnv_sv/sawfish_discover/{sample}_{type}/expected.copy.number.bed"
            if config.get("sawfish_discover", {}).get("expected_cn", False)
            else []
        ),
        gc_mpack="cnv_sv/sawfish_discover/{sample}_{type}/genome.gclevels.mpack",
        max_dp="cnv_sv/sawfish_discover/{sample}_{type}/max.depth.bed",
        stats_json="cnv_sv/sawfish_discover/{sample}_{type}/run.stats.json",
        gcbias_mpack="cnv_sv/sawfish_discover/{sample}_{type}/sample.gcbias.mpack",
    output:
        stats_json=temp("cnv_sv/sawfish_joint_call/{sample}_{type}/run_stats.json"),
        log=temp("cnv_sv/sawfish_joint_call/{sample}_{type}/sawfish.log"),
        sample_vcf=temp("cnv_sv/sawfish_joint_call/{sample}_{type}.vcf.gz"),
        sample_tbi=temp("cnv_sv/sawfish_joint_call/{sample}_{type}.vcf.gz.tbi"),
        gt_vcf=temp("cnv_sv/sawfish_joint_call/{sample}_{type}/genotyped.sv.vcf.gz"),
        gt_tbi=temp("cnv_sv/sawfish_joint_call/{sample}_{type}/genotyped.sv.vcf.gz.tbi"),
        sr_json=(
            temp("cnv_sv/sawfish_joint_call/{sample}_{type}/supporting_reads.json.gz")
            if config.get("sawfish_joint_call", {}).get("supporting_reads", False)
            else []
        ),
    params:
        extra=config.get("sawfish_joint_call", {}).get("extra", ""),
        in_dir="cnv_sv/sawfish_discover/{sample}_{type}/",
        out_dir="cnv_sv/sawfish_joint_call/{sample}_{type}",
        suporting_reads=(
            f"--report-supporting-reads" if config.get("sawfish_joint_call", {}).get("supporting_reads", False) else ""
        ),
    log:
        "cnv_sv/sawfish_joint_call/{sample}_{type}.joint_call.log",
    benchmark:
        repeat(
            "cnv_sv/sawfish_joint_call/{sample}_{type}.joint_call.benchmark.tsv",
            config.get("sawfish_joint_call", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("sawfish_joint_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sawfish_joint_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sawfish_joint_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sawfish_joint_call", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sawfish_joint_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sawfish_joint_call", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sawfish_joint_call", {}).get("container", config["default_container"])
    message:
        "{rule}: Run sawish joint call on the ouput of sawfish discover"
    shell:
        "sawfish joint-call "
        "{params.extra} "
        "--threads {threads} "
        "--sample {params.in_dir} "
        "--output-dir {params.out_dir} "
        "--report-supporting-reads &> {log} && "
        "cp {output.gt_vcf} {output.sample_vcf} && "
        "cp {output.gt_tbi} {output.sample_tbi} "
