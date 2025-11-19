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
        repeat("cnv_sv/melt/{sample}_{type}.melt.benchmark.tsv", config.get("melt", {}).get("benchmark_repeats", 1))
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
