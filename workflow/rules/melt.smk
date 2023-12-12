__author__ = "Padraic Corcoran and Magdalena"
__copyright__ = "Copyright 2023, Padraic Corcoran and Magdalena"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule melt:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam", #  /D23-05906_TC83_230818/D23-05906_TC83_230818-dedup.bam 
        ref=config.get("reference", {}).get("fasta", ""),
        bed_file=config.get("melt", {}).get("include_bed", ""),
        mei_file=config.get("melt", {}).get("include_mei", ""),
    output:
        alu="cnv_sv/melt/Results_single_{sample}_{type}/ALU.final_comp.vcf",
        line1="cnv_sv/melt/Results_single_{sample}_{type}/LINE1.final_comp.vcf",
        sva="cnv_sv/melt/Results_single_{sample}_{type}/SVA.final_comp.vcf",
    params:
        folder="cnv_sv/melt/Results_single_{sample}_{type}",  # /melt_res/Results_single_D23-05906_TC83_230818
        extra=config.get("melt", {}).get("extra", ""),
    log:
        "cnv_sv/melt/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/melt/{sample}_{type}.output.benchmark.tsv",
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
        "{rule}: Do stuff on cnv_sv/{rule}/{wildcards.sample}_{wildcards.type}.input"
    shell:
        "(/melt.sh '-Xmx{resources.mem_mb}m' Single " 
        "-w {params.folder} "
        "-bamfile {input.bam} " 
        "-h {input.ref} "
        "-t {input.mei_file} "
        "-n {input.bed_file} "
        "{params.extra}) "
        "&> {log}"       


