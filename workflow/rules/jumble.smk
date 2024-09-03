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
        cnvkit_segments=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.cns"),
        cnvkit_bins=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.cnr"),
        counts=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.counts.RDS"),
        jumble_segments=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.segments.csv"),
        jumble_bins=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.jumble.RDS"),
        output_dir=temp(directory("cnv_sv/jumble_run/{sample}_{type}/")),
        png=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.png"),
        snps=temp("cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.jumble.snps.RDS"),
    params:
        reference=config.get("jumble", {}).get("normal_reference", ""),
    log:
        "cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.output.benchmark.tsv",
            config.get("jumble_run", {}).get("benchmark_repeats", 1)
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
