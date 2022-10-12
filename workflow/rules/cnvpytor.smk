__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2022, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvpytor_readdepth:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        vcf="parabricks/pbrun_deepvariant/{sample}.vcf",
    output:
        pytor=temp("cnv_sv/cnvpytor/{sample}.pytor"),
    params:
        extra=config.get("cnvpytor", {}).get("extra", ""),
    log:
        "cnv_sv/cnvpytor/{sample}_rd.log"
    benchmark:
        repeat("cnv_sv/cnvpytor/{sample}_rd.output.benchmark.tsv",
        config.get("cnvpytor", {}).get("benchmark_repeats", 1))
    threads: config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvpytor", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvpytor", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvpytor", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvpytor", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvpytor", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvpytor.yaml"
    message:
       "{rule}: Run cnvpytor calls based on read depth for {wildcards.sample}_{wildcards.type}"
    shell:
        """cnvpytor -root {output.pytor} -rd {input.bam} &&
        cnvpytor -root {output.pytor} -his 1000 10000 100000 &&
        cnvpytor -root {output.pytor} -partition 1000 10000 100000 &&
        cnvpytor -root {output.pytor} -call 1000 10000 100000 &&
        cnvpytor -root {output.pytor} -snp {input.vcf} -sample {wildcards.sample}_N &&
        cnvpytor -root {output.pytor} -mask_snps &&
        cnvpytor -root {output.pytor} -baf 10000 100000 &> {log}"""

# Joint segmentation & caller
# cnvpytor -root file.pytor -call combined 10000

rule cnvpytor_filter:
    input:
        pytor=temp("cnv_sv/cnvpytor/{sample}.pytor")
    output:
        xls="cnv_sv/cnvpytor/{sample}.xls",
        vcf="cnv_sv/cnvpytor/{sample}.vcf",
    params:
        extra=config.get("cnvpytor", {}).get("extra", ""),
    log:
        "cnv_sv/cnvpytor/{sample}_filter.log"
    benchmark:
        repeat("cnv_sv/cnvpytor/{sample}_filter.output.benchmark.tsv",
        config.get("cnvpytor", {}).get("benchmark_repeats", 1))
    threads: config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvpytor", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvpytor", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvpytor", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvpytor", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvpytor", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvpytor.yaml"
    message:
       "{rule}: Filter cnvpytor calls for {wildcards.sample}"
    shell:
        """cnvpytor -root {input.pytor} -view 100000  &&
        print calls &&
        set Q0_range 0 0.5 &&
        set size_range 100000 inf &&
        print calls &&
        set p_range 0 0.00001 &&
        set print_filename {output.xls} &&
        print calls &&
        set print_filename {output.vcf} &&
        print calls &&
        &> {log}"""

rule cnvpytor_annotate:
    input:
        xls="cnv_sv/cnvpytor/{sample}.xls",
        vcf="cnv_sv/cnvpytor/{sample}.vcf",
    output:
        tsv="cnv_sv/cnvpytor/{sample}.tsv",
    params:
        extra=config.get("cnvpytor", {}).get("extra", ""),
    log:
        "cnv_sv/cnvpytor/{sample}_annotate.log"
    benchmark:
        repeat("cnv_sv/cnvpytor/{sample}_annotate.output.benchmark.tsv",
        config.get("cnvpytor", {}).get("benchmark_repeats", 1))
    threads: config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvpytor", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvpytor", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvpytor", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvpytor", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvpytor", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvpytor.yaml"
    message:
       "{rule}: Annotate cnvpytor calls for {wildcards.sample}"
    shell:
        """cnvpytor -root {input.pytor} -view 100000 &&
        set Q0_range 0 0.5 &&
        set size_range 100000 inf &&
        set print_filename {output.tsv} &&
        set annotate &&
        print calls &&
        &> {log}"""


rule cnvpytor_plot:
    input:
        tsv="cnv_sv/cnvpytor/{sample}.tsv",
        pytor=temp("cnv_sv/cnvpytor/{sample}.pytor")
    output:
        tsv="cnv_sv/cnvpytor/{sample}.png",
    params:
        extra=config.get("cnvpytor", {}).get("extra", ""),
    log:
        "cnv_sv/cnvpytor/{sample}_plot.log"
    benchmark:
        repeat("cnv_sv/cnvpytor/{sample}_plot.output.benchmark.tsv",
        config.get("cnvpytor", {}).get("benchmark_repeats", 1))
    threads: config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvpytor", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvpytor", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvpytor", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvpytor", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvpytor", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvpytor.yaml"
    message:
       "{rule}: Plot cnvpytor calls for {wildcards.sample}"
    shell:
        """echo "rdstat" | cnvpytor -root {input.pytor} -view 100000 -o {output.png} &&
        cnvpytor -root file.pytor -view 100000 <<ENDL \
        set rd_use_mask \
        set markersize 1 \
        set grid vertical \
        set output_filename {output.png} \
        manhattan \
        circular \
        rd \
        likelihood \
        baf \
        ENDL &&
        cnvpytor -root {input.pytor} -view 100000 < script.spytor &&
        &> {log}"""
