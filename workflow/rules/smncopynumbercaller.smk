__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule smn_caller:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
    output:
        json=temp("cnv_sv/smn_caller/{sample}_{type}.json"),
        tsv=temp("cnv_sv/smn_caller/{sample}_{type}.tsv"),
    params:
        extra=config.get("smn_caller", {}).get("extra", ""),
        genome=config.get("smn_caller", {}).get("genome_version", ""),
        outdir=lambda wildcards, output: "{}".format(os.path.split(output.tsv)[0]),
        prefix=lambda wildcards, output: "{}_{}".format(wildcards.sample, wildcards.type),
    log:
        "cnv_sv/smn_caller/{sample}_{type}.tsv.log",
    benchmark:
        repeat(
            "cnv_sv/smn_caller/{sample}_{type}.tsv.benchmark.tsv",
            config.get("smn_caller", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("smn_caller", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("smn_caller", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("smn_caller", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("smn_caller", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("smn_caller", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("smn_caller", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("smn_caller", {}).get("container", config["default_container"])
    conda:
        "../envs/smncopynumbercaller.yaml"
    message:
        "{rule}: Call SMN1 and SMN2 copynumber on {input.bam} using SMNCopyNumberCaller"
    shell:
        "find $(pwd)/{input.bam} > {input.bam}.manifest.txt && "
        "smn_caller.py --manifest {input.bam}.manifest.txt "
        "--genome {params.genome} "
        "--prefix {params.prefix} "     
        "--outDir {params.outdir} "      
        "--threads {resources.threads} &> {log}"      


rule smn_charts:
    input:
        json="cnv_sv/smn_caller/{sample}_{type}.json"
    output:
        outdir=temp(directory("cnv_sv/smn_charts/{sample}_{type}")),
        pdf=temp("cnv_sv/smn_charts/{sample}_{type}.pdf"),
    params:
        extra=config.get("smn_charts", {}).get("extra", ""),
        genome=config.get("smn_charts", {}).get("genome_version", ""),
        outdir=lambda wildcards, output: "{}".format(os.path.split(output.pdf)[0]),
    log:
        "cnv_sv/smn_charts/{sample}_{type}.pdf.log",
    benchmark:
        repeat(
            "cnv_sv/smn_charts/{sample}_{type}.pdf.benchmark.tsv",
            config.get("smn_charts", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("smn_charts", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("smn_charts", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("smn_charts", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("smn_charts", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("smn_charts", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("smn_charts", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("smn_charts", {}).get("container", config["default_container"])
    conda:
        "../envs/smncopynumbercaller.yaml"
    message:
        "{rule}: Visualisation of SMNCopyNumberCaller result in {input.json}"
    shell:
        "python $(whereis smn_charts.py | awk '{{print $2}}') -s {input.json} " 
        "-o {params.outdir} &> {log}"        