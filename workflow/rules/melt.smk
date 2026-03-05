__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2025, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule melt:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        bed=config.get("melt", {}).get("bed", ""),
        mei=config.get("melt", {}).get("mei", ""),
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        alu=temp("cnv_sv/melt/{sample}_{type}.ALU.final_comp.vcf"),
        hervk=temp("cnv_sv/melt/{sample}_{type}.HERVK.final_comp.vcf"),
        line1=temp("cnv_sv/melt/{sample}_{type}.LINE1.final_comp.vcf"),
        sva=temp("cnv_sv/melt/{sample}_{type}.SVA.final_comp.vcf"),
        tmpdir=temp(directory("cnv_sv/melt/{sample}_{type}")),
    params:
        extra=config.get("melt", {}).get("extra", ""),
    log:
        "cnv_sv/melt/{sample}_{type}.melt.output.log",
    benchmark:
        repeat(
            "cnv_sv/melt/{sample}_{type}.melt.benchmark.tsv",
            config.get("melt", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("melt", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("melt", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("melt", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("melt", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("melt", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("melt", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("melt", {}).get("container", config["default_container"])
    message:
        "{rule}: Calling mobile elements in {input.bam} with MELT"
    shell:
        """
        export JAVA_OPTS="-Xmx{resources.mem_mb}m"
        (MELT \
        Single \
        -bamfile {input.bam} \
        -h {input.ref} \
        -t {input.mei} \
        -n {input.bed} \
        -w {output.tmpdir} \
        {params.extra} && \
        cp {output.tmpdir}/ALU.final_comp.vcf {output.alu} && \
        cp {output.tmpdir}/HERVK.final_comp.vcf {output.hervk} && \
        cp {output.tmpdir}/LINE1.final_comp.vcf {output.line1} && \
        cp {output.tmpdir}/SVA.final_comp.vcf {output.sva}) \
        &> {log}
        """


rule melt_concat:
    input:
        alu_vcf="cnv_sv/melt/{sample}_{type}.ALU.final_comp.vcf.gz",
        alu_tbi="cnv_sv/melt/{sample}_{type}.ALU.final_comp.vcf.gz.tbi",
        hervk_vcf="cnv_sv/melt/{sample}_{type}.HERVK.final_comp.vcf.gz",
        hervk_tbi="cnv_sv/melt/{sample}_{type}.HERVK.final_comp.vcf.gz.tbi",
        line1_vcf="cnv_sv/melt/{sample}_{type}.LINE1.final_comp.vcf.gz",
        line1_tbi="cnv_sv/melt/{sample}_{type}.LINE1.final_comp.vcf.gz.tbi",
        sva_vcf="cnv_sv/melt/{sample}_{type}.SVA.final_comp.vcf.gz",
        sva_tbi="cnv_sv/melt/{sample}_{type}.SVA.final_comp.vcf.gz.tbi",
    output:
        vcf="cnv_sv/melt/{sample}_{type}.concat.vcf.gz",
    log:
        "cnv_sv/melt/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/melt/{sample}_{type}.vcf.benchmark.tsv",
            config.get("melt_concat", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("melt_concat", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("melt_concat", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("melt_concat", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("melt_concat", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("melt_concat", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("melt_concat", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("melt_concat", {}).get("container", config["default_container"])
    message:
        "{rule}: Concatenating MELT variant files for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        bcftools concat \\
        -a {input.alu_vcf} {input.hervk_vcf} {input.line1_vcf} {input.sva_vcf} \\
        | bcftools sort -Oz -o {output.vcf} \\
        &> {log}
        """


rule melt_vcf:
    input:
        vcf="cnv_sv/melt/{sample}_{type}.concat.vcf.gz",
    output:
        vcf="cnv_sv/melt/{sample}_{type}.vcf",
    params:
        extra=config.get("melt_vcf", {}).get("extra", ""),
    log:
        "cnv_sv/melt/{sample}_{type}.melt_vcf.log",
    benchmark:
        repeat("cnv_sv/melt/{sample}_{type}.vcf.benchmark.tsv", config.get("melt_vcf", {}).get("benchmark_repeats", 1))
    threads: config.get("melt_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("melt_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("melt_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("melt_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("melt_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("melt_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("melt_vcf", {}).get("container", config["default_container"])
    message:
        "{rule}: normalise the VCF file {input.vcf}"
    script:
        "../scripts/melt_vcf.py"
