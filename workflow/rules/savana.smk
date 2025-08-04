__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2025, Andrei Guliaev"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"


rule savana_pb_to:
    input:
        bam="annotation/whatshap_haplotag/{sample}_{type}.haplotagged.bam",
        bai="annotation/whatshap_haplotag/{sample}_{type}.haplotagged.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
        fai=config.get("reference", {}).get("fai", ""),
    output:
        dir=temp(directory("cnv_sv/savana_pb_to/{sample}_{type}")),
    params:
        min_support=config.get("savana_pb_to", {}).get("min_support", 10),
        extra=config.get("savana_pb_to", {}).get("extra", ""),
    log:
        "cnv_sv/savana_pb_to/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/savana_pb_to/{sample}_{type}.output.benchmark.tsv", config.get("savana_pb_to", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("savana_pb_to", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("savana_pb_to", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("savana_pb_to", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("savana_pb_to", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("savana_pb_to", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("savana_pb_to", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("savana_pb_to", {}).get("container", config["default_container"])
    message:
        "{rule}: run savana on {input.bam}"
    shell:
        "savana to --tumour {input.bam} "
        "--outdir {output.dir} "
        "--ref {input.ref} "
        "--ref_index {input.fai} "
        "--pb --threads {threads} "
        "--min_support {params.min_support} "
        "{params.extra} &> {log}"


rule savana_pb_tn:
    input:
        bam_n="annotation/whatshap_haplotag/{sample}_N.haplotagged.bam",
        bai_n="annotation/whatshap_haplotag/{sample}_N.haplotagged.bam",
        bam_t="annotation/whatshap_haplotag/{sample}_T.haplotagged.bam",
        bai_t="annotation/whatshap_haplotag/{sample}_T.haplotagged.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
        fai=config.get("reference", {}).get("fai", ""),
    output:
        dir=temp(directory("cnv_sv/savana_pb_tn/{sample}")),
    params:
        min_support=config.get("savana_pb_tn", {}).get("min_support", 10),
        extra=config.get("savana_pb_tn", {}).get("extra", ""),
    log:
        "cnv_sv/savana_pb_tn/{sample}.output.log",
    benchmark:
        repeat("cnv_sv/savana_pb_tn/{sample}.output.benchmark.tsv", config.get("savana_pb_tn", {}).get("benchmark_repeats", 1))
    threads: config.get("savana_pb_tn", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("savana_pb_tn", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("savana_pb_tn", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("savana_pb_tn", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("savana_pb_tn", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("savana_pb_tn", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("savana_pb_tn", {}).get("container", config["default_container"])
    message:
        "{rule}: run savana on {input.bam_n} and {input.bam_t}"
    shell:
        "savana --tumour {input.bam_t} "
        "--normal {input.bam_n} "
        "--outdir {output.dir} "
        "--ref {input.ref} "
        "--ref_index {input.fai} "
        "--pb --threads {threads} "
        "--min_support {params.min_support} "
        "{params.extra} &> {log}"


rule savana_ont_to:
    input:
        bam="annotation/whatshap_haplotag/{sample}_{type}.haplotagged.bam",
        bai="annotation/whatshap_haplotag/{sample}_{type}.haplotagged.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
        fai=config.get("reference", {}).get("fai", ""),
    output:
        dir=temp(directory("cnv_sv/savana_ont_to/{sample}_{type}")),
    params:
        extra=config.get("savana_ont_to", {}).get("extra", ""),
    log:
        "cnv_sv/savana_ont_to/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/savana_ont_to/{sample}_{type}.output.benchmark.tsv",
            config.get("savana_ont_to", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("savana_ont_to", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("savana_ont_to", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("savana_ont_to", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("savana_ont_to", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("savana_ont_to", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("savana_ont_to", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("savana_ont_to", {}).get("container", config["default_container"])
    message:
        "{rule}: run savana on {input.bam}"
    shell:
        "savana to --tumour {input.bam} "
        "--outdir {output.dir} "
        "--ref {input.ref} "
        "--ref_index {input.fai} "
        "--ont --threads {threads} "
        "{params.extra} &> {log}"


rule savana_ont_tn:
    input:
        bam_n="annotation/whatshap_haplotag/{sample}_N.haplotagged.bam",
        bai_n="annotation/whatshap_haplotag/{sample}_N.haplotagged.bam",
        bam_t="annotation/whatshap_haplotag/{sample}_T.haplotagged.bam",
        bai_t="annotation/whatshap_haplotag/{sample}_T.haplotagged.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
        fai=config.get("reference", {}).get("fai", ""),
    output:
        dir=temp(directory("cnv_sv/savana_ont_tn/{sample}")),
    params:
        extra=config.get("savana_ont_tn", {}).get("extra", ""),
    log:
        "cnv_sv/savana_ont_tn/{sample}.output.log",
    benchmark:
        repeat("cnv_sv/savana_ont_tn/{sample}.output.benchmark.tsv", config.get("savana_ont_tn", {}).get("benchmark_repeats", 1))
    threads: config.get("savana_ont_tn", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("savana_ont_tn", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("savana_ont_tn", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("savana_ont_tn", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("savana_ont_tn", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("savana_ont_tn", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("savana_ont_tn", {}).get("container", config["default_container"])
    message:
        "{rule}: run savana on {input.bam_n} and {input.bam_t}"
    shell:
        "savana --tumour {input.bam_t} "
        "--normal {input.bam_n} "
        "--outdir {output.dir} "
        "--ref {input.ref} "
        "--ref_index {input.fai} "
        "--ont --threads {threads} "
        "{params.extra} &> {log}"
