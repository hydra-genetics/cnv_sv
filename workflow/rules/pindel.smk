__author__ = "Martin Rippin"
__copyright__ = "Copyright 2022, Martin Rippin"
__email__ = "martin.rippin@scilifelab.uu.se"
__license__ = "GPL-3"


rule pindel_generate_config:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        metrics="qc/picard_collect_multiple_metrics/{sample}_{type}.insert_size_metrics",
    output:
        config=temp("cnv_sv/pindel/{sample}_{type}.cfg"),
    log:
        "cnv_sv/pindel/{sample}_{type}_pindel_generate_config.log",
    benchmark:
        repeat(
            "cnv_sv/pindel/{sample}_{type}_pindel_generate_config.benchmark.tsv",
            config.get("pindel_generate_config", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("pindel_generate_config", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("pindel_generate_config", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pindel_generate_config", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pindel_generate_config", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("pindel_generate_config", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pindel_generate_config", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("pindel_generate_config", {}).get("container", config["default_container"])
    message:
        "{rule}: produce config for {wildcards.sample} {wildcards.type}"
    script:
        "../scripts/generate_pindel_config.py"


rule pindel_call:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        config="cnv_sv/pindel/{sample}_{type}.cfg",
        ref=config["reference"]["fasta"],
        include_bed=config.get("pindel_call", {}).get("include_bed", []),
        exclude_bed=config.get("pindel_call", {}).get("exclude_bed", []),
    output:
        pindel=temp(
            expand(
                "cnv_sv/pindel/{{sample}}_{{type}}_{ext}",
                ext=[
                    "BP",
                    "CloseEndMapped",
                    "D",
                    "INT_final",
                    "INV",
                    "LI",
                    "RP",
                    "SI",
                    "TD",
                ],
            )
        ),
    params:
        extra=config.get("pindel_call", {}).get("extra", ""),
    log:
        "cnv_sv/pindel/{sample}_{type}_pindel_call.log",
    benchmark:
        repeat(
            "cnv_sv/pindel/{sample}_{type}_pindel_call.benchmark.tsv",
            config.get("pindel_call", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("pindel_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("pindel_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pindel_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pindel_call", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("pindel_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pindel_call", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("pindel_call", {}).get("container", config["default_container"])
    message:
        "{rule}: detect breakpoints in {wildcards.sample} {wildcards.type}"
    wrapper:
        "v2.6.0/bio/pindel/call"


rule pindel2vcf:
    input:
        pindel=expand(
            "cnv_sv/pindel/{{sample}}_{{type}}_{ext}",
            ext=[
                "BP",
                "CloseEndMapped",
                "D",
                "INT_final",
                "INV",
                "LI",
                "RP",
                "SI",
                "TD",
            ],
        ),
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("cnv_sv/pindel_vcf/{sample}_{type}.no_contig.vcf"),
    params:
        extra=config.get("pindel2vcf", {}).get("extra", ""),
        refdate=config.get("pindel2vcf", {}).get("refdate", "20131217"),
        refname=config.get("pindel2vcf", {}).get("refname", "'Genome Reference Consortium Human Build 38'"),
    log:
        "cnv_sv/pindel_vcf/{sample}_{type}.no_contig.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/pindel/{sample}_{type}.no_contig.vcf.benchmark.tsv",
            config.get("pindel2vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("pindel2vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("pindel2vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pindel2vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pindel2vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("pindel2vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pindel2vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("pindel2vcf", {}).get("container", config["default_container"])
    message:
        "{rule}: convert pindel output to vcf for {wildcards.sample}_{wildcards.type}.no_contig"
    wrapper:
        "v2.6.0/bio/pindel/pindel2vcf"


rule pindel_update_vcf:
    input:
        fasta=config["reference"]["fasta"],
        vcf="cnv_sv/pindel_vcf/{sample}_{type}.no_contig.vcf",
    output:
        vcf=temp("cnv_sv/pindel_vcf/{sample}_{type}.no_tc.vcf"),
        samplename=temp("cnv_sv/pindel_vcf/{sample}_{type}.samplename.txt"),
    params:
        extra=config.get("pindel_update_vcf", {}).get("extra", ""),
    log:
        "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/pindel_vcf/{sample}_{type}.no_tc.vcf.benchmark.tsv",
            config.get("pindel_update_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("pindel_update_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("pindel_update_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pindel_update_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pindel_update_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("pindel_update_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pindel_update_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("pindel_update_vcf", {}).get("container", config["default_container"])
    message:
        "{rule}: update cnv_sv/pindel/{wildcards.sample}_{wildcards.type}.no_contig.vcf to include contigs and correct samplename"
    shell:
        "(echo -e '{wildcards.sample}\n' > {output.samplename} && "
        "picard UpdateVcfSequenceDictionary "
        "-INPUT {input.vcf} "
        "-SD {input.fasta} "
        "-OUTPUT /dev/stdout | bcftools reheader -s {output.samplename} -o {output.vcf} - )&> {log}"
