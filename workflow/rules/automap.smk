__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2023, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu."
__license__ = "GPL-3"


rule automap:
    input:
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vcf",
    output:
        pdf=temp("cnv_sv/automap/{sample}_{type}/{sample}_{type}.HomRegions.pdf"),
        tsv=temp("cnv_sv/automap/{sample}_{type}/{sample}_{type}.HomRegions.tsv"),
    params:
        extra=config.get("automap", {}).get("extra", ""),
        build=config.get("automap", {}).get("build", ""),
        dir=temp(directory(lambda w, output: os.path.dirname(os.path.dirname(output[0])))),
    log:
        "cnv_sv/automap/{sample}_{type}.output.log",
    benchmark:
        repeat("cnv_sv/automap/{sample}_{type}.output.benchmark.tsv", config.get("automap", {}).get("benchmark_repeats", 1))
    threads: config.get("automap", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("automap", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("automap", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("automap", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("automap", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("automap", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("automap", {}).get("container", config["default_container"])
    conda:
        "../envs/automap.yaml"
    message:
        "{rule}: Finding ROH regions {output.tsv}"
    shell:
        "automap "
        "--vcf {input.vcf} "
        "--out {params.dir} "
        "--genome {params.build} "
        "{params.extra} &> {log}"
