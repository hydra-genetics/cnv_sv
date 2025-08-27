__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2024, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule sawfish_discover:
    input:
        bam=lambda wildcards: get_input_bam(wildcards)[0],
        bai=lambda wildcards: get_input_bam(wildcards)[1],
        maf=(
            "snv_indels/deepvariant/{sample}_{type}.merged.vcf.gz"
            if config.get("sawfish_discover", {}).get("maf", False)
            else []
        ),
        maf_tbi=(
            "snv_indels/deepvariant/{sample}_{type}.merged.vcf.gz.tbi"
            if config.get("sawfish_discover", {}).get("maf", False)
            else []
        ),
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        asm_bed=temp("cnv_sv/sawfish_discover/{sample}_{type}/assembly.regions.bed"),
        bcf=temp("cnv_sv/sawfish_discover/{sample}_{type}/candidate.sv.bcf"),
        csi=temp("cnv_sv/sawfish_discover/{sample}_{type}/candidate.sv.bcf.csi"),
        contig_bam=temp("cnv_sv/sawfish_discover/{sample}_{type}/contig.alignment.bam"),
        contig_csi=temp("cnv_sv/sawfish_discover/{sample}_{type}/contig.alignment.bam.csi"),
        cluster_bed=temp("cnv_sv/sawfish_discover/{sample}_{type}/debug.breakpoint_clusters.bed"),
        cn_bed=(
            "cnv_sv/sawfish_discover/{sample}_{type}/expected.copy.number.bed"
            if config.get("sawfish_discover", {}).get("expected_cn", False)
            else []
        ),
        cn_bdg=(
            "cnv_sv/sawfish_discover/{sample}_{type}/copynum.bedgraph"
            if config.get("sawfish_discover", {}).get("disable_cnv", False)
            else []
        ),
        cn_mpack=(
            "cnv_sv/sawfish_discover/{sample}_{type}/copynum.mpack"
            if config.get("sawfish_discover", {}).get("disable_cnv", False)
            else []
        ),
        dp_mpack=temp("cnv_sv/sawfish_discover/{sample}_{type}/depth.mpack"),
        gc_mpack=(
            "cnv_sv/sawfish_discover/{sample}_{type}/genome.gclevels.mpack"
            if config.get("sawfish_discover", {}).get("disable_cnv", False)
            else []
        ),
        gcbias_mpack=(
            "cnv_sv/sawfish_discover/{sample}_{type}/sample.gcbias.mpack"
            if config.get("sawfish_discover", {}).get("disable_cnv", False)
            else []
        ),
        max_dp=temp("cnv_sv/sawfish_discover/{sample}_{type}/max.depth.bed"),
        refine_txt=temp("cnv_sv/sawfish_discover/{sample}_{type}/debug.cluster.refinement.txt"),
        settings_json=temp("cnv_sv/sawfish_discover/{sample}_{type}/discover.settings.json"),
        stats_json=temp("cnv_sv/sawfish_discover/{sample}_{type}/run.stats.json"),
    params:
        extra=config.get("sawfish_discover", {}).get("extra", ""),
        expected_cn=get_expected_cn,
        disable_cnv=("--disable-cnv" if config.get("sawfish_discover", {}).get("disable_cnv", False) else ""),
        out_dir="cnv_sv/sawfish_discover/{sample}_{type}",
        maf=(f"--maf {input.maf} " if config.get("sawfish_discover", {}).get("maf", False) else ""),
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
        "{params.maf} "
        "--threads {threads} "
        "--ref {input.ref} "
        "--bam {input.bam} "
        "--output-dir {params.out_dir} "
        "{params.disable_cnv} "
        "{params.expected_cn} &> {log}"


rule sawfish_joint_call_single:
    input:
        bam=lambda wildcards: get_input_bam(wildcards)[0],
        bai=lambda wildcards: get_input_bam(wildcards)[1],
        ref=config.get("reference", {}).get("fasta", ""),
        asm_bed="cnv_sv/sawfish_discover/{sample}_{type}/assembly.regions.bed",
        bcf="cnv_sv/sawfish_discover/{sample}_{type}/candidate.sv.bcf",
        csi="cnv_sv/sawfish_discover/{sample}_{type}/candidate.sv.bcf.csi",
        contig_bam="cnv_sv/sawfish_discover/{sample}_{type}/contig.alignment.bam",
        contig_csi="cnv_sv/sawfish_discover/{sample}_{type}/contig.alignment.bam.csi",
        cluster_bed="cnv_sv/sawfish_discover/{sample}_{type}/debug.breakpoint_clusters.bed",
        refine_txt="cnv_sv/sawfish_discover/{sample}_{type}/debug.cluster.refinement.txt",
        dp_mpack="cnv_sv/sawfish_discover/{sample}_{type}/depth.mpack",
        settings_json="cnv_sv/sawfish_discover/{sample}_{type}/discover.settings.json",
        cn_bed=(
            "cnv_sv/sawfish_discover/{sample}_{type}/expected.copy.number.bed"
            if config.get("sawfish_discover", {}).get("expected_cn", False)
            else []
        ),
        cn_bdg=(
            "cnv_sv/sawfish_discover/{sample}_{type}/copynum.bedgraph"
            if config.get("sawfish_discover", {}).get("disable_cnv", False)
            else []
        ),
        cn_mpack=(
            "cnv_sv/sawfish_discover/{sample}_{type}/copynum.mpack"
            if config.get("sawfish_discover", {}).get("disable_cnv", False)
            else []
        ),
        gc_mpack=(
            "cnv_sv/sawfish_discover/{sample}_{type}/genome.gclevels.mpack"
            if config.get("sawfish_discover", {}).get("disable_cnv", False)
            else []
        ),
        gcbias_mpack=(
            "cnv_sv/sawfish_discover/{sample}_{type}/sample.gcbias.mpack"
            if config.get("sawfish_discover", {}).get("disable_cnv", False)
            else []
        ),
        max_dp="cnv_sv/sawfish_discover/{sample}_{type}/max.depth.bed",
        stats_json="cnv_sv/sawfish_discover/{sample}_{type}/run.stats.json",
    output:
        cn_bdg=(
            temp("cnv_sv/sawfish_joint_call/{sample}_{type}/samples/sample0001_{sample}_{type}/copynum.bedgraph")
            if config.get("sawfish_joint_call", {}).get("disable_cnv", False)
            else []
        ),
        dp_bw=temp("cnv_sv/sawfish_joint_call/{sample}_{type}/samples/sample0001_{sample}_{type}/depth.bw"),
        gcbias_bw=(
            temp("cnv_sv/sawfish_joint_call/{sample}_{type}/samples/sample0001_{sample}_{type}/gc_bias_corrected_depth.bw")
            if config.get("sawfish_joint_call", {}).get("disable_cnv", False)
            else []
        ),
        gt_vcf=temp("cnv_sv/sawfish_joint_call/{sample}_{type}/genotyped.sv.vcf.gz"),
        gt_tbi=temp("cnv_sv/sawfish_joint_call/{sample}_{type}/genotyped.sv.vcf.gz.tbi"),
        log=temp("cnv_sv/sawfish_joint_call/{sample}_{type}/sawfish.log"),
        maf=(
            temp("cnv_sv/sawfish_joint_call/{sample}_{type}/samples/sample_0001_{sample}_{type}/maf.bw")
            if config.get("sawfish_joint_call", {}).get("maf", False)
            else []
        ),
        stats_json=temp("cnv_sv/sawfish_joint_call/{sample}_{type}/run.stats.json"),
        sample_vcf=temp("cnv_sv/sawfish_joint_call/{sample}_{type}.vcf.gz"),
        sample_tbi=temp("cnv_sv/sawfish_joint_call/{sample}_{type}.vcf.gz.tbi"),
        sr_json=(
            temp("cnv_sv/sawfish_joint_call/{sample}_{type}/supporting_reads.json.gz")
            if config.get("sawfish_joint_call", {}).get("supporting_reads", False)
            else []
        ),
    params:
        extra=config.get("sawfish_joint_call", {}).get("extra", ""),
        in_dir=lambda w, input: os.path.dirname(input.settings_json),
        out_dir=lambda w, output: output[0].rsplit("/", 3)[0],
        supporting_reads=(
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
        "{rule}: Run sawfish joint call on the single sample ouput of sawfish discover"
    shell:
        "sawfish joint-call "
        "{params.extra} "
        "--threads {threads} "
        "--sample {params.in_dir} "
        "--output-dir {params.out_dir} "
        "{params.supporting_reads} &> {log} && "
        "cp {output.gt_vcf} {output.sample_vcf} && "
        "cp {output.gt_tbi} {output.sample_tbi} "
