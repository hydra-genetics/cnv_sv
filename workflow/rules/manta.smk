__author__ = "Martin Rippin, Jonas Almlöf, Jessika Nordin"
__copyright__ = "Copyright 2021, Martin Rippin, Jonas Almlöf"
__email__ = "martin.rippin@scilifelab.uu.se, jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule manta_config_tn:
    input:
        bam_t="alignment/samtools_merge_bam/{sample}_T.bam",
        bai_t="alignment/samtools_merge_bam/{sample}_T.bam.bai",
        bam_n="alignment/samtools_merge_bam/{sample}_N.bam",
        bai_n="alignment/samtools_merge_bam/{sample}_N.bam.bai",
        ref=config["reference"]["fasta"],
    output:
        scrpt=temp("cnv_sv/manta_run_workflow_tn/{sample}/runWorkflow.py"),
    params:
        extra=config.get("manta_config_tn", {}).get("extra", ""),
    log:
        "cnv_sv/manta_config_tn/{sample}/runWorkflow.py.log",
    benchmark:
        repeat(
            "cnv_sv/manta_config_tn/{sample}/runWorkflow.py.benchmark.tsv",
            config.get("manta_config_tn", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("manta_config_tn", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("manta_config_tn", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("manta_config_tn", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("manta_config_tn", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("manta_config_tn", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("manta_config_tn", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("manta_config_tn", {}).get("container", config["default_container"])
    message:
        "{rule}: generate manta runWorkflow.py for {wildcards.sample}"
    shell:
        "configManta.py "
        "--tumorBam={input.bam_t} "
        "--normalBam={input.bam_n} "
        "--referenceFasta={input.ref} "
        "{params.extra} "
        "--runDir=cnv_sv/manta_run_workflow_tn/{wildcards.sample} &> {log}"


rule manta_run_workflow_tn:
    input:
        bam_t="alignment/samtools_merge_bam/{sample}_T.bam",
        bai_t="alignment/samtools_merge_bam/{sample}_T.bam.bai",
        bam_n="alignment/samtools_merge_bam/{sample}_N.bam",
        bai_n="alignment/samtools_merge_bam/{sample}_N.bam.bai",
        ref=config["reference"]["fasta"],
        scrpt="cnv_sv/manta_run_workflow_tn/{sample}/runWorkflow.py",
    output:
        cand_si_vcf=temp("cnv_sv/manta_run_workflow_tn/{sample}/results/variants/candidateSmallIndels.vcf.gz"),
        cand_si_tbi=temp("cnv_sv/manta_run_workflow_tn/{sample}/results/variants/candidateSmallIndels.vcf.gz.tbi"),
        cand_sv_vcf=temp("cnv_sv/manta_run_workflow_tn/{sample}/results/variants/candidateSV.vcf.gz"),
        cand_sv_tbi=temp("cnv_sv/manta_run_workflow_tn/{sample}/results/variants/candidateSV.vcf.gz.tbi"),
        dipl_sv_vcf=temp("cnv_sv/manta_run_workflow_tn/{sample}/results/variants/diploidSV.vcf.gz"),
        dipl_sv_tbi=temp("cnv_sv/manta_run_workflow_tn/{sample}/results/variants/diploidSV.vcf.gz.tbi"),
        som_sv_vcf=temp("cnv_sv/manta_run_workflow_tn/{sample}/results/variants/somaticSV.vcf.gz"),
        som_sv_tbi=temp("cnv_sv/manta_run_workflow_tn/{sample}/results/variants/somaticSV.vcf.gz.tbi"),
        wrk_dir=temp(directory("cnv_sv/manta_run_workflow_tn/{sample}/workspace")),
    log:
        "cnv_sv/manta_run_workflow_tn/{sample}/manta_tn.log",
    benchmark:
        repeat(
            "cnv_sv/manta_run_workflow_tn/{sample}/manta_tn.benchmark.tsv",
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
    message:
        "{rule}: use manta to call sv in {wildcards.sample}"
    shell:
        "{input.scrpt} "
        "-j {threads} "
        "-g unlimited &> {log}"


rule manta_config_t:
    input:
        bam_t="alignment/samtools_merge_bam/{sample}_T.bam",
        bai_t="alignment/samtools_merge_bam/{sample}_T.bam.bai",
        ref=config["reference"]["fasta"],
    output:
        scrpt=temp("cnv_sv/manta_run_workflow_t/{sample}/runWorkflow.py"),
    params:
        extra=config.get("manta_config_t", {}).get("extra", ""),
    log:
        "cnv_sv/manta_config_t/{sample}/runWorkflow.py.log",
    benchmark:
        repeat(
            "cnv_sv/manta_config_t/{sample}/runWorkflow.py.benchmark.tsv",
            config.get("manta_config_t", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("manta_config_t", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("manta_config_t", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("manta_config_t", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("manta_config_t", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("manta_config_t", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("manta_config_t", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("manta_config_t", {}).get("container", config["default_container"])
    message:
        "{rule}: generate manta runWorkflow.py for {wildcards.sample}"
    shell:
        "configManta.py "
        "--tumorBam={input.bam_t} "
        "--referenceFasta={input.ref} "
        "{params.extra} "
        "--runDir=cnv_sv/manta_run_workflow_t/{wildcards.sample} &> {log}"


rule manta_run_workflow_t:
    input:
        bam_t="alignment/samtools_merge_bam/{sample}_T.bam",
        bai_t="alignment/samtools_merge_bam/{sample}_T.bam.bai",
        ref=config["reference"]["fasta"],
        scrpt="cnv_sv/manta_run_workflow_t/{sample}/runWorkflow.py",
    output:
        cand_si_vcf=temp("cnv_sv/manta_run_workflow_t/{sample}/results/variants/candidateSmallIndels.vcf.gz"),
        cand_si_tbi=temp("cnv_sv/manta_run_workflow_t/{sample}/results/variants/candidateSmallIndels.vcf.gz.tbi"),
        cand_sv_vcf=temp("cnv_sv/manta_run_workflow_t/{sample}/results/variants/candidateSV.vcf.gz"),
        cand_sv_tbi=temp("cnv_sv/manta_run_workflow_t/{sample}/results/variants/candidateSV.vcf.gz.tbi"),
        tum_sv_vcf=temp("cnv_sv/manta_run_workflow_t/{sample}/results/variants/tumorSV.vcf.gz"),
        tum_sv_tbi=temp("cnv_sv/manta_run_workflow_t/{sample}/results/variants/tumorSV.vcf.gz.tbi"),
        wrk_dir=temp(directory("cnv_sv/manta_run_workflow_t/{sample}/workspace")),
    log:
        "cnv_sv/manta_run_workflow_t/{sample}/manta_t.log",
    benchmark:
        repeat(
            "cnv_sv/manta_run_workflow_t/{sample}/manta_t.benchmark.tsv",
            config.get("manta_run_workflow_t", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("manta_run_workflow_t", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("manta_run_workflow_t", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("manta_run_workflow_t", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("manta_run_workflow_t", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("manta_run_workflow_t", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("manta_run_workflow_t", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("manta_run_workflow_t", {}).get("container", config["default_container"])
    message:
        "{rule}: use manta to call sv in {wildcards.sample}"
    shell:
        "{input.scrpt} "
        "-j {threads} "
        "-g unlimited &> {log}"


rule manta_config_n:
    input:
        bam_n="alignment/samtools_merge_bam/{sample}_N.bam",
        bai_n="alignment/samtools_merge_bam/{sample}_N.bam.bai",
        ref=config["reference"]["fasta"],
    output:
        scrpt=temp("cnv_sv/manta_run_workflow_n/{sample}/runWorkflow.py"),
    params:
        extra=config.get("manta_config_n", {}).get("extra", ""),
    log:
        "cnv_sv/manta_config_n/{sample}/runWorkflow.py.log",
    benchmark:
        repeat(
            "cnv_sv/manta_config_n/{sample}/runWorkflow.py.benchmark.tsv",
            config.get("manta_config_n", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("manta_config_n", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("manta_config_n", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("manta_config_n", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("manta_config_n", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("manta_config_n", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("manta_config_n", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("manta_config_n", {}).get("container", config["default_container"])
    message:
        "{rule}: generate manta runWorkflow.py for {wildcards.sample}"
    shell:
        "configManta.py "
        "--bam={input.bam_n} "
        "--referenceFasta={input.ref} "
        "{params.extra} "
        "--runDir=cnv_sv/manta_run_workflow_n/{wildcards.sample} &> {log}"


rule manta_run_workflow_n:
    input:
        bam_n="alignment/samtools_merge_bam/{sample}_N.bam",
        bai_n="alignment/samtools_merge_bam/{sample}_N.bam.bai",
        ref=config["reference"]["fasta"],
        scrpt="cnv_sv/manta_run_workflow_n/{sample}/runWorkflow.py",
    output:
        cand_si_vcf=temp("cnv_sv/manta_run_workflow_n/{sample}/results/variants/candidateSmallIndels.vcf.gz"),
        cand_si_tbi=temp("cnv_sv/manta_run_workflow_n/{sample}/results/variants/candidateSmallIndels.vcf.gz.tbi"),
        cand_sv_vcf=temp("cnv_sv/manta_run_workflow_n/{sample}/results/variants/candidateSV.vcf.gz"),
        cand_sv_tbi=temp("cnv_sv/manta_run_workflow_n/{sample}/results/variants/candidateSV.vcf.gz.tbi"),
        wrk_dir=temp(directory("cnv_sv/manta_run_workflow_n/{sample}/workspace")),
    log:
        "cnv_sv/manta_run_workflow_n/{sample}/manta_run_workflow_n.output.log",
    benchmark:
        repeat(
            "cnv_sv/manta_run_workflow_n/{sample}/manta_run_workflow_n.output.benchmark.tsv",
            config.get("manta_run_workflow_n", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("manta_run_workflow_n", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("manta_run_workflow_n", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("manta_run_workflow_n", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("manta_run_workflow_n", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("manta_run_workflow_n", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("manta_run_workflow_n", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("manta_run_workflow_n", {}).get("container", config["default_container"])
    message:
        "{rule}: use manta to call sv in {wildcards.sample}"
    shell:
        "{input.scrpt} "
        "-j {threads} "
        "-g unlimited &> {log}"
