__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule gatk_collect_read_counts:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        interval=config.get("reference", {}).get("design_intervals_gatk_cnv", ""),
    output:
        temp("cnv_sv/gatk_collect_read_counts/{sample}_{type}.counts.hdf5"),
    params:
        extra=config.get("gatk_collect_read_counts", {}).get("extra", ""),
        mergingRule="OVERLAPPING_ONLY",
    log:
        "cnv_sv/gatk_collect_read_counts/{sample}_{type}.counts.hdf5.log",
    benchmark:
        repeat(
            "cnv_sv/gatk_collect_read_counts/{sample}_{type}.counts.hdf5.benchmark.tsv",
            config.get("gatk_collect_read_counts", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_collect_read_counts", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_collect_read_counts", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_collect_read_counts", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_collect_read_counts", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_collect_read_counts", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_collect_read_counts", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_collect_read_counts", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: use gatk_cnv to obtain cnv_sv/gatk_collect_read_counts/{wildcards.sample}_{wildcards.type}.counts.hdf5"
    shell:
        "(gatk --java-options '-Xmx4g' CollectReadCounts "
        "-I {input.bam} "
        "-L {input.interval} "
        "--interval-merging-rule {params.mergingRule} "
        "{params.extra} "
        "-O {output}) &> {log}"


rule gatk_collect_allelic_counts:
    input:
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        interval=config.get("gatk_collect_allelic_counts", {}).get("SNP_interval", ""),
        ref=config["reference"]["fasta"],
    output:
        temp("cnv_sv/gatk_collect_allelic_counts/{sample}_{type}.clean.allelicCounts.tsv"),
    params:
        extra=config.get("gatk_collect_allelic_counts", {}).get("extra", ""),
    log:
        "cnv_sv/gatk_collect_allelic_counts/{sample}_{type}.clean.allelicCounts.tsv.log",
    benchmark:
        repeat(
            "cnv_sv/gatk_collect_allelic_counts/{sample}_{type}.clean.allelicCounts.tsv.benchmark.tsv",
            config.get("gatk_collect_allelic_counts", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_collect_allelic_counts", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_collect_allelic_counts", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_collect_allelic_counts", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_collect_allelic_counts", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_collect_allelic_counts", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_collect_allelic_counts", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_collect_allelic_counts", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: use gatk_cnv to obtain cnv_sv/gatk_collect_allelic_counts/{wildcards.sample}_{wildcards.type}.clean.allelicCounts.tsv"
    shell:
        "(gatk --java-options '-Xmx4g' CollectAllelicCounts "
        "-L {input.interval} "
        "-I {input.bam} "
        "-R {input.ref} "
        "-O {output}"
        "{params.extra}) &> {log}"


rule gatk_denoise_read_counts:
    input:
        hdf5PoN=config.get("gatk_denoise_read_counts", {}).get("normal_reference", ""),
        hdf5Tumor="cnv_sv/gatk_collect_read_counts/{sample}_{type}.counts.hdf5",
    output:
        denoisedCopyRatio=temp("cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv"),
        stdCopyRatio=temp("cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.standardizedCR.tsv"),
    params:
        extra=config.get("gatk_denoise_read_counts", {}).get("extra", ""),
    log:
        "cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv.log",
    benchmark:
        repeat(
            "cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv.benchmark.tsv",
            config.get("gatk_denoise_read_counts", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_denoise_read_counts", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_denoise_read_counts", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_denoise_read_counts", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_denoise_read_counts", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_denoise_read_counts", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_denoise_read_counts", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_denoise_read_counts", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: use gatk_cnv to obtain cnv_sv/gatk_denoise_read_counts/{wildcards.sample}_{wildcards.type}.clean.denoisedCR.tsv"
    shell:
        "(gatk --java-options '-Xmx4g' DenoiseReadCounts "
        "-I {input.hdf5Tumor} "
        "--count-panel-of-normals {input.hdf5PoN} "
        "--standardized-copy-ratios {output.stdCopyRatio} "
        "--denoised-copy-ratios {output.denoisedCopyRatio} "
        "{params.extra}) &> {log}"


rule gatk_model_segments:
    input:
        allelicCounts="cnv_sv/gatk_collect_allelic_counts/{sample}_{type}.clean.allelicCounts.tsv",
        denoisedCopyRatio="cnv_sv/gatk_denoise_read_counts/{sample}_{type}.clean.denoisedCR.tsv",
    output:
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelFinal.seg"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.cr.seg"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.af.igv.seg"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.cr.igv.seg"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.hets.tsv"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelBegin.cr.param"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelBegin.af.param"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelBegin.seg"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelFinal.af.param"),
        temp("cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelFinal.cr.param"),
    params:
        extra=config.get("gatk_model_segments", {}).get("extra", ""),
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        outprefix="{sample}_{type}.clean",
    log:
        "cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelFinal.seg.log",
    benchmark:
        repeat(
            "cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelFinal.seg.benchmark.tsv",
            config.get("gatk_model_segments", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_model_segments", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_model_segments", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_model_segments", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_model_segments", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_model_segments", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_model_segments", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_model_segments", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: use gatk to obtain cnv_sv/gatk_model_segments/{wildcards.sample}_{wildcards.type}.clean.modelFinal.seg"
    shell:
        "(gatk --java-options '-Xmx4g' ModelSegments "
        "--denoised-copy-ratios {input.denoisedCopyRatio} "
        "--allelic-counts {input.allelicCounts} "
        "--output {params.outdir} "
        "--output-prefix {params.outprefix}"
        "{params.extra}) &> {log}"


rule gatk_call_copy_ratio_segments:
    input:
        "cnv_sv/gatk_model_segments/{sample}_{type}.clean.cr.seg",
    output:
        segments=temp("cnv_sv/gatk_call_copy_ratio_segments/{sample}_{type}.clean.calledCNVs.seg"),
        igv_segments=temp("cnv_sv/gatk_call_copy_ratio_segments/{sample}_{type}.clean.calledCNVs.igv.seg"),
    params:
        extra=config.get("gatk_call_copy_ratio_segments", {}).get("extra", ""),
    log:
        "cnv_sv/gatk_call_copy_ratio_segments/{sample}_{type}.clean.calledCNVs.seg.log",
    benchmark:
        repeat(
            "cnv_sv/gatk_call_copy_ratio_segments/{sample}_{type}.clean.calledCNVs.seg.benchmark.tsv",
            config.get("gatk_call_copy_ratio_segments", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_call_copy_ratio_segments", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_call_copy_ratio_segments", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_call_copy_ratio_segments", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_call_copy_ratio_segments", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_call_copy_ratio_segments", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_call_copy_ratio_segments", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_call_copy_ratio_segments", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: use gatk_to obtain cnv_sv/gatk_call_copy_ratio_segments/{wildcards.sample}_{wildcards.type}.clean.calledCNVs.seg"
    shell:
        "(gatk --java-options '-Xmx4g' CallCopyRatioSegments "
        "--input {input} "
        "--output {output.segments} "
        "{params.extra}) &> {log}"


rule gatk_cnv_to_vcf:
    input:
        segment="cnv_sv/gatk_model_segments/{sample}_{type}.clean.modelFinal.seg",
    output:
        vcf=temp("cnv_sv/gatk_vcf/{sample}_{type}.vcf"),
    params:
        dup_limit=config.get("gatk_cnv_vcf", {}).get("dup_limit", 2.5),
        het_del_limit=config.get("gatk_cnv_vcf", {}).get("het_del_limit", 1.5),
        hom_del_limit=config.get("gatk_cnv_vcf", {}).get("hom_del_limit", 0.5),
        sample_id="{sample}_{type}",
        tc=lambda wildcards: get_sample(samples, wildcards)["tumor_content"],
    log:
        "cnv_sv/gatk_cnv_vcf/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/gatk_cnv_vcf/{sample}_{type}.vcf.benchmark.tsv",
            config.get("gatk_cnv_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_cnv_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_cnv_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_cnv_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_cnv_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_cnv_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_cnv_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_cnv_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/gatk.yaml"
    message:
        "{rule}: export gatk cnv segments into vcf in cnv_sv/gatk_cnv_vcf/{wildcards.sample}_{wildcards.type}.vcf"
    script:
        "../scripts/gatk_cnv_to_vcf.py"
