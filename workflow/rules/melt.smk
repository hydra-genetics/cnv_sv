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
        alu=temp("cnv_sv/melt/{sample}_{type}/ALU.final_comp.vcf.gz"),
        line1=temp("cnv_sv/melt/{sample}_{type}/LINE1.final_comp.vcf.gz"),
        sva=temp("cnv_sv/melt/{sample}_{type}/SVA.final_comp.vcf.gz"),
        tmpdir=temp(directory("cnv_sv/melt/{sample}_{type}")),
    params:
        extra=config.get("melt", {}).get("extra", ""),
    log:
        "cnv_sv/melt/{sample}_{type}/melt.output.log",
    benchmark:
        repeat(
            "cnv_sv/melt/{sample}_{type}/melt.output.benchmark.tsv",
            config.get("melt", {}).get("benchmark_repeats", 1)
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
        (java -Xmx{resources.mem_mb}m -jar \
        /projects/wp3/Software/MELTv2.2.2/MELT.jar \
        Single \
        -bamfile {input.bam} \
        -h {input.ref} \
        -t {input.mei} \
        -n {input.bed} \
        -w {output.tmpdir} \
        {params.extra}) \
        &> {log}
        """
