__author__ = "Padraic Corcoran, Jessika Nordin"
__copyright__ = "Copyright 2023, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule upd_regions:
    input:
        vcf="snv_indels/glnexus/{sample}_{type}.vep_annotated.vcf.gz",
    output:
        regions="cnv_sv/upd/{sample}_{type}.upd_regions.bed",
    params:
        father=lambda wildcards: get_parent_samples(wildcards, "father"),
        mother=lambda wildcards: get_parent_samples(wildcards, "mother"),
        proband="{sample}_{type}",
        extra=config.get("upd", {}).get("extra", ""),
    log:
        "cnv_sv/upd/{sample}_{type}.upd_regions.bed.log",
    benchmark:
        repeat(
            "cnv_sv/upd/{sample}_{type}.upd_regions.bed.benchmark.tsv",
            config.get("upd", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("upd", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("upd", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("upd", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("upd", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("upd", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("upd", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("upd", {}).get("container", config["default_container"])
    message:
        "{rule}: Use upd on a trio vcf {input.vcf} to find regions of uniparental disomy"
    shell:
        "(upd "
        "--vcf {input.vcf} "
        "--proband {params.proband} "
        "--mother {params.mother} "
        "--father {params.father} "
        "{params.extra} "
        "regions "
        "--out {output.regions}) "
        "&> {log}"


rule upd_sites:
    input:
        vcf="snv_indels/glnexus/{sample}_{type}.vep_annotated.vcf.gz",
    output:
        sites="cnv_sv/upd/{sample}_{type}.upd_sites.bed",
    params:
        father=lambda wildcards: get_parent_samples(wildcards, "father"),
        mother=lambda wildcards: get_parent_samples(wildcards, "mother"),
        proband="{sample}_{type}",
        extra=config.get("upd", {}).get("extra", ""),
    log:
        "cnv_sv/upd/{sample}_{type}.upd_sites.bed.log",
    benchmark:
        repeat(
            "cnv_sv/upd/{sample}_{type}.upd_sites.bed.benchmark.tsv",
            config.get("upd", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("upd", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("upd", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("upd", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("upd", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("upd", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("upd", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("upd", {}).get("container", config["default_container"])
    message:
        "{rule}: Use upd on a trio vcf {input.vcf} to find uniparental disomy informative sites"
    shell:
        "(upd "
        "--vcf {input.vcf} "
        "--proband {params.proband} "
        "--mother {params.mother} "
        "--father {params.father} "
        "{params.extra} "
        "sites "
        "--out {output.sites}) "
        "&> {log}"
