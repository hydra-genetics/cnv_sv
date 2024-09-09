__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2024, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule jumble_run:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.germline.vcf",
    output:
        cnvkit_segments=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.bam.cns"),
        cnvkit_bins=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.bam.cnr"),
        counts=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.bam.counts.RDS"),
        dna_copy_segments=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.bam_dnacopy.seg"),
        jumble_segments=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.bam.segments.csv"),
        jumble_bins=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.bam.jumble.RDS"),
        output_dir=temp(directory("cnv_sv/jumble_run/{sample}_{type}/")),
        png=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.bam.png"),
        snps=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.bam.jumble.snps.RDS"),
    params:
        reference=config.get("jumble_run", {}).get("normal_reference", ""),
    log:
        "cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.output.benchmark.tsv",
            config.get("jumble_run", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("jumble_run", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("jumble_run", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("jumble_run", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("jumble_run", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("jumble_run", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("jumble_run", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("jumble_run", {}).get("container", config["default_container"])
    message:
        "{rule}: Call CNVs with jumble on {input.bam}"
    shell:
        "(Rscript /Jumble/jumble-run.R "
        "-r {params.reference} "
        "-b {input.bam} "
        "-v {input.vcf} "
        "-o {output.output_dir} &> {log}"


rule jumble_cnvkit_call:
    input:
        segments="cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.bam.cns",
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.germline.vcf",
        tc_file=get_tc_file,
    output:
        segment=temp("cnv_sv/jumble_cnvkit_call/{sample}_{type}.{tc_method}.loh.cns"),
    params:
        extra=config.get("jumble_cnvkit_call", {}).get("extra", ""),
        purity=get_tc,
    log:
        "cnv_sv/jumble_cnvkit_call/{sample}_{type}.{tc_method}.loh.cns.log",
    benchmark:
        repeat(
            "cnv_sv/jumble_cnvkit_call/{sample}_{type}.{tc_method}.loh.cns.benchmark.tsv",
            config.get("jumble_cnvkit_call", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("jumble_cnvkit_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("jumble_cnvkit_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("jumble_cnvkit_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("jumble_cnvkit_call", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("jumble_cnvkit_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("jumble_cnvkit_call", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("jumble_cnvkit_call", {}).get("container", config["default_container"])
    message:
        "{rule}: call cnvs with loh info into {output.segment}"
    wrapper:
        "v3.3.6/bio/cnvkit/call"


rule jumble_vcf:
    input:
        segment="cnv_sv/jumble_cnvkit_call/{sample}_{type}.{tc_method}.loh.cns",
        tc_file=get_tc_file,
    output:
        vcf=temp("cnv_sv/jumble_vcf/{sample}_{type}.{tc_method}.vcf"),
    params:
        dup_limit=config.get("jumble_vcf", {}).get("dup_limit", 2.5),
        het_del_limit=config.get("jumble_vcf", {}).get("het_del_limit", 1.5),
        hom_del_limit=config.get("jumble_vcf", {}).get("hom_del_limit", 0.5),
        sample_id="{sample}_{type}",
        tc=get_tc,
    log:
        "cnv_sv/jumble_vcf/{sample}_{type}.{tc_method}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/jumble_vcf/{sample}_{type}.{tc_method}.vcf.output.benchmark.tsv",
            config.get("jumble_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("jumble_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("jumble_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("jumble_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("jumble_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("jumble_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("jumble_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("jumble_vcf", {}).get("container", config["default_container"])
    message:
        "{rule}: export cnvkit segments into vcf in {output.vcf}"
    script:
        "../scripts/cnvkit_vcf.py"
