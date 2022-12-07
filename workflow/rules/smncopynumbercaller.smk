__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule smn_manifest:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
    output:
        manifest=temp("cnv_sv/smn_caller/{sample}_{type}_manifest.txt")
    params:
        extra=config.get("smn_caller", {}).get("extra", ""),
    log:
        "cnv_sv/smn_caller/{sample}_{type}_manifest.txt.log",
    benchmark:
        repeat(
            "cnv_sv/smn_caller/{sample}_{type}_manifest.txt.benchmark.tsv",
            config.get("smn_caller", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("smn_manifest", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("smn_manifest", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("smn_manifest", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("smn_manifest", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("smn_manifest", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("smn_manifest", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("smn_manifest", {}).get("container", config["default_container"])
    conda:
        "../envs/smncopynumbercaller.yaml"
    message:
        "{rule}: Generate the manifest file for SMNCopyNumberCaller"
    script:
        "../scripts/smncopynumbercaller_manifest.py"      


rule smn_caller:
    input:
       "cnv_sv/smn_caller/{sample}_{type}_manifest.txt",
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
        "{rule}: Call SMN1 and SMN2 copynumber on {input} using SMNCopyNumberCaller"
    shell:
        "smn_caller.py --manifest {input} "
        "--genome {params.genome} "
        "--prefix {params.prefix} "     
        "--outDir {params.outdir} "      
        "--threads {resources.threads} &> {log}"      


rule smn_charts:
    input:
        json="cnv_sv/smn_caller/{sample}_{type}.json"
    output:
        pdf=temp("cnv_sv/smn_charts/smn_{sample}_{type}.pdf"),
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