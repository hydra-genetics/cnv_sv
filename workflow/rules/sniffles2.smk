__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2024, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule sniffles2_call:
    input:
        bam=lambda wildcards: get_input_aligned_bam(wildcards, config)[0],
        bai=lambda wildcards: get_input_aligned_bam(wildcards, config)[1],
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf=temp("cnv_sv/sniffles2_call/{sample}_{type}.vcf"),
        snf=temp("cnv_sv/sniffles2_call/{sample}_{type}.snf"),
    params:
        sample_id=lambda wildcards, output: "{}_{}".format(wildcards.sample, wildcards.type),
        tandem_repeats=get_tr_bed,
        extra=config.get("sniffles2_call", {}).get("extra", ""),
    log:
        "cnv_sv/sniffles2_call/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/sniffles2_call/{sample}_{type}.output.benchmark.tsv",
            config.get("sniffles2_call", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("sniffles2_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sniffles2_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sniffles2_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sniffles2_call", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sniffles2_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sniffles2_call", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sniffles2_call", {}).get("container", config["default_container"])
    message:
        "{rule}: Calls SVs on {input.bam} with sniffles"
    shell:
        "sniffles -i {input.bam} "
        "--reference {input.ref} "
        "-t {threads} "
        "--sample-id {params.sample_id} "
        "{params.tandem_repeats} "
        "{params.extra} "
        "--vcf {output.vcf} "
        "--snf {output.snf} &> {log}"


rule sniffles2_joint_call:
    input:
        snfs=lambda wildcards: [
            f"cnv_sv/sniffles2_call/{wildcards.sample}_{t}.snf"
            for t in get_unit_types(units, wildcards.sample)
        ],
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf=temp("cnv_sv/sniffles2_joint_call/{sample}.vcf.gz"),
    params:
        extra=config.get("sniffles2_joint_call", {}).get("extra", ""),
    log:
        "cnv_sv/sniffles2_joint_call/{sample}.vcf.gz.log",
    benchmark:
        repeat(
            "cnv_sv/sniffles2_joint_call/{sample}.vcf.gz.benchmark.tsv",
            config.get("sniffles2_joint_call", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("sniffles2_joint_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sniffles2_joint_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sniffles2_joint_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sniffles2_joint_call", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sniffles2_joint_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sniffles2_joint_call", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sniffles2_joint_call", {}).get("container", config["default_container"])
    message:
        "{rule}: Sniffles2 joint call for {wildcards.sample} ({input.snfs})"
    shell:
        """
        (sniffles \
            --input {input.snfs} \
            --vcf {output.vcf} \
            --reference {input.ref} \
            --threads {threads} \
            {params.extra}) 2> {log}
        """
