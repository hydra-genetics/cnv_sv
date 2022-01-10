__author__ = "Martin Rippin"
__copyright__ = "Copyright 2021, Martin Rippin"
__email__ = "martin.rippin@scilifelab.uu.se"
__license__ = "GPL-3"


rule config_manta_tn:
    input:
        bam_t="alignment/merge_bam/{sample}_T.bam",
        bai_t="alignment/merge_bam/{sample}_T.bam.bai",
        bam_n="alignment/merge_bam/{sample}_N.bam",
        bai_n="alignment/merge_bam/{sample}_N.bam.bai",
        ref=config["reference"]["fasta"],
    output:
        scrpt=temp("cnv_sv/manta/{sample}/runWorkflow.py"),
    log:
        "cnv_sv/manta/{sample}/run_workflow.py.log",
    benchmark:
        repeat(
            "cnv_sv/manta/{sample}/run_workflow.py.benchmark.tsv",
            config.get("config_manta_tn", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("config_manta_tn", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("config_manta_tn", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("config_manta_tn", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("config_manta_tn", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("config_manta_tn", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("config_manta_tn", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("config_manta_tn", {}).get("container", config["default_container"])
    conda:
        "../envs/config_manta_tn.yaml"
    message:
        "{rule}: Generate manta runWorkflow.py for {wildcards.sample}"
    shell:
        "configManta.py "
        "--tumorBam={input.bam_t} "
        "--normalBam={input.bam_n} "
        "--referenceFasta={input.ref} "
        "--runDir=cnv_sv/manta &> {log}"
