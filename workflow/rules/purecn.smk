__author__ = "Erik Demitz-Helin"
__copyright__ = "Copyright 2022, Erik Demitz-Helin"
__email__ = "erik.demitz-helin@gu.se"
__license__ = "GPL-3"


rule purecn_coverage:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
    output:
        purecn=temp(
            expand(
                "cnv_sv/purecn_coverage/{{sample}}_{{type}}{ext}",
                ext=[
                    "_coverage.txt.gz",
                    "_coverage_loess.txt.gz",
                    "_coverage_loess.png",
                    "_coverage_loess_qc.txt",
                ],
            )
        ),
    params:
        intervals=config.get("purecn_coverage", {}).get("intervals", ""),
        extra=config.get("purecn_coverage", {}).get("extra", ""),
    log:
        "cnv_sv/purecn_coverage/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/purecn_coverage/{sample}_{type}.output.benchmark.tsv",
            config.get("purecn_coverage", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("purecn_coverage", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("purecn_coverage", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("purecn_coverage", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("purecn_coverage", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("purecn_coverage", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("purecn_coverage", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("purecn_coverage", {}).get("container", config["default_container"])
    message:
        "{rule}: calculate coverage for {wildcards.sample}"
    shell:
        "(Rscript $PURECN/Coverage.R "
        "--out-dir=cnv_sv/purecn_coverage "
        "--bam={input.bam} "
        "--intervals={params.intervals} "
        "--cores={threads} "
        "{params.extra}) &> {log}"


checkpoint purecn:
    input:
        unpack(get_purecn_inputs),
        vcf="snv_indels/gatk_mutect2/{sample}_{type}.merged.unfiltered.bcftools_annotated.vcf.gz",
        tbi="snv_indels/gatk_mutect2/{sample}_{type}.merged.unfiltered.bcftools_annotated.vcf.gz.tbi",
    output:
        csv=temp("cnv_sv/purecn/temp/{sample}_{type}/{sample}_{type}.csv"),
        outdir=temp(directory("cnv_sv/purecn/temp/{sample}_{type}/")),
    params:
        fun_segmentation=config.get("purecn", {}).get("fun_segmentation", ""),
        genome=config.get("purecn", {}).get("genome", ""),
        interval_padding=config.get("purecn", {}).get("interval_padding", ""),
        extra=get_purecn_extra,
    log:
        "cnv_sv/purecn/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/purecn/{sample}_{type}.output.benchmark.tsv",
            config.get("purecn", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("purecn", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("purecn", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("purecn", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("purecn", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("purecn", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("purecn", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("purecn", {}).get("container", config["default_container"])
    message:
        "{rule}: Quantify purity/ploidy for {wildcards.sample}_{wildcards.type}"
    shell:
        "(Rscript $PURECN/PureCN.R "
        "--tumor={input.tumor} "
        "--vcf={input.vcf} "
        "--genome={params.genome} "
        "--fun-segmentation={params.fun_segmentation} "
        "--interval-padding={params.interval_padding} "
        "--sampleid={wildcards.sample}_{wildcards.type} "
        "--out={output.outdir} "
        "--force --seed=1337 "
        "{params.extra} || touch {output.csv}) &> {log}"


rule purecn_copy_output:
    input:
        file_dir=lambda wildcards: checkpoints.purecn.get(**wildcards).output.outdir,
        files=lambda wildcards: f"{checkpoints.purecn.get(** wildcards).output.outdir}/{{sample}}_{{type}}{{suffix}}",
    output:
        files=temp("cnv_sv/purecn/{sample}_{type}{suffix}"),
    wildcard_constraints:
        suffix="|".join(
            f"({s})"
            for s in [
                ".csv",
                ".rds",
                ".pdf",
                "_dnacopy.seg",
                "_chromosomes.pdf",
                "_segmentation.pdf",
                "_local_optima.pdf",
                "_variants.csv",
                "_loh.csv",
            ]
        ),
    log:
        "cnv_sv/purecn/{sample}_{type}{suffix}.output.log",
    benchmark:
        repeat(
            "cnv_sv/purecn/{sample}_{type}{suffix}.purecn_copy_output.benchmark.tsv",
            config.get("purecn_copy_output", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("purecn_copy_output", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("purecn_copy_output", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("purecn_copy_output", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("purecn_copy_output", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("purecn_copy_output", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("purecn_copy_output", {}).get("time", config["default_resources"]["time"]),
    container:
        config["default_container"]
    message:
        "{rule}: Copy {input.files} to {output.files}"
    shell:
        "cp {input.files} {output.files}"


rule purecn_purity_file:
    input:
        csv="cnv_sv/purecn/{sample}_{type}.csv",
    output:
        purity=temp("cnv_sv/purecn_purity_file/{sample}_{type}.purity.txt"),
    params:
        min_purity=config.get("purecn_purity_file", {}).get("min_purity", "0.2"),
        missing_purity_value=config.get("purecn_purity_file", {}).get("missing_purity_value", "1"),
    log:
        "cnv_sv/purecn_purity_file/{sample}_{type}.purity.txt.log",
    benchmark:
        repeat(
            "cnv_sv/purecn_purity_file/{sample}_{type}.purity.txt.benchmark.tsv",
            config.get("purecn_purity_file", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("purecn_purity_file", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("purecn_purity_file", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("purecn_purity_file", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("purecn_purity_file", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("purecn_purity_file", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("purecn_purity_file", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("purecn_purity_file", {}).get("container", config["default_container"])
    message:
        "{rule}: Extract purity value for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        (awk -F',' 'FNR==2 {{ print ($2 > {params.min_purity} ? $2 : {params.min_purity}) }} \
        END{{if (NR==1) {{ print {params.missing_purity_value} }} }}' \
        {input.csv} > {output.purity}) &> {log}
        """
