__author__ = "Martin Rippin"
__copyright__ = "Copyright 2021, Martin Rippin"
__email__ = "martin.rippin@scilifelab.uu.se"
__license__ = "GPL-3"


rule manta_run_workflow_tn:
    input:
        bam_t="alignment/merge_bam/{sample}_T.bam",
        bai_t="alignment/merge_bam/{sample}_T.bam.bai",
        bam_n="alignment/merge_bam/{sample}_N.bam",
        bai_n="alignment/merge_bam/{sample}_N.bam.bai",
        ref=config["reference"]["fasta"],
        scrpt="cnv_sv/manta/{sample}/runWorkflow.py",
    output:
        cand_si_vcf=temp("cnv_sv/manta/{sample}/results/variants/candidateSmallIndels.vcf.gz"),
        cand_si_tbi=temp("cnv_sv/manta/{sample}/results/variants/candidateSmallIndels.vcf.gz.tbi"),
        cand_sv_vcf=temp("cnv_sv/manta/{sample}/results/variants/candidateSV.vcf.gz"),
        cand_sv_tbi=temp("cnv_sv/manta/{sample}/results/variants/candidateSV.vcf.gz.tbi"),
        dipl_sv_vcf=temp("cnv_sv/manta/{sample}/results/variants/diploidSV.vcf.gz"),
        dipl_sv_tbi=temp("cnv_sv/manta/{sample}/results/variants/diploidSV.vcf.gz.tbi"),
        som_sv_vcf="cnv_sv/manta/{sample}/results/variants/somaticSV.vcf.gz",
        som_sv_tbi="cnv_sv/manta/{sample}/results/variants/somaticSV.vcf.gz.tbi",
        wrk_dir=temp(directory("cnv_sv/manta/{sample}/workspace")),
    log:
        "cnv_sv/manta/{sample}/manta_tn.log",
    benchmark:
        repeat(
            "cnv_sv/manta/{sample}/manta_tn.benchmark.tsv",
            config.get("manta_run_workflow_tn", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("manta_run_workflow_tn", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("manta_run_workflow_tn", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("manta_run_workflow_tn", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("manta_run_workflow_tn", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("manta_run_workflow_tn", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("manta_run_workflow_tn", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("manta_run_workflow_tn", {}).get("container", config["default_container"])
    conda:
        "../envs/manta_run_workflow_tn.yaml"
    message:
        "{rule}: Use manta to call sv in {wildcards.sample}"
    shell:
        "{input.scrpt} "
        "-j {threads} "
        "-g unlimited &> {log}"
