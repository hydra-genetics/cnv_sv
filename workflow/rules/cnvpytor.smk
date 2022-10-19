__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2022, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvpytor_readdepth:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vcf.gz",
    output:
        pytor=temp("cnv_sv/cnvpytor/{sample}_{type}.pytor"),
    params:
        extra=config.get("cnvpytor", {}).get("extra", ""),
    log:
        "cnv_sv/cnvpytor/{sample}_{type}_rd.log"
    benchmark:
        repeat("cnv_sv/cnvpytor/{sample}_{type}_rd.output.benchmark.tsv",
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
        cnvpytor -root {output.pytor} -snp {input.vcf} -sample {wildcards.sample}_{wildcards.type} &&
        cnvpytor -root {output.pytor} -mask_snps &&
        cnvpytor -root {output.pytor} -baf 1000 10000 100000 &&
        cnvpytor -root {output.pytor} -call combined 1000 10000 100000
        &> {log}"""

# Joint segmentation & caller
# cnvpytor -root file.pytor -call combined 10000

rule cnvpytor_filter:
    input:
        pytor="cnv_sv/cnvpytor/{sample}_{type}.pytor",
    output:
        vcf="cnv_sv/cnvpytor/{sample}_{type}.vcf",
        filtvcf="cnv_sv/cnvpytor/{sample}_{type}_filtered.vcf",
    params:
        extra=config.get("cnvpytor", {}).get("extra", ""),
    log:
        "cnv_sv/cnvpytor/{sample}_{type}_filter.log"
    benchmark:
        repeat("cnv_sv/cnvpytor/{sample}_{type}_filter.output.benchmark.tsv",
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
       "{rule}: Filter cnvpytor calls for {wildcards.sample}_{wildcards.type}"
    shell:
        """singularity run cnvpytor_latest.sif cnvpytor -root {input.pytor} -view 1000 <<-ENDL
        set print_filename {output.vcf}
        print calls
        set Q0_range 0 0.5
        set p_range 0 0.01
        set pN_range 0 0.05
        set dG_range 100000 inf
        set print_filename {output.filtvcf}
        print calls
        ENDL
        &> {log}"""


## Might be added later, code not tested ##
#rule cnvpytor_plot:
#    input:
#        tsv="cnv_sv/cnvpytor/{sample}.tsv",
#        pytor="cnv_sv/cnvpytor/{sample}_{type}.pytor",
#    output:
#        tsv="cnv_sv/cnvpytor/{sample}.png",
#    params:
#        extra=config.get("cnvpytor", {}).get("extra", ""),
#    log:
#        "cnv_sv/cnvpytor/{sample}_plot.log"
#    benchmark:
#        repeat("cnv_sv/cnvpytor/{sample}_plot.output.benchmark.tsv",
#        config.get("cnvpytor", {}).get("benchmark_repeats", 1))
#    threads: config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"])
#    resources:
#        mem_mb=config.get("cnvpytor", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
#        mem_per_cpu=config.get("cnvpytor", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
#        partition=config.get("cnvpytor", {}).get("partition", config["default_resources"]["partition"]),
#        threads=config.get("cnvpytor", {}).get("threads", config["default_resources"]["threads"]),
#        time=config.get("cnvpytor", {}).get("time", config["default_resources"]["time"]),
#    container:
#        config.get("cnvpytor", {}).get("container", config["default_container"])
#    conda:
#        "../envs/cnvpytor.yaml"
#    message:
##    shell:
#        """echo "rdstat" | cnvpytor -root {input.pytor} -view 100000 -o {output.png} &&
#        cnvpytor -root file.pytor -view 100000 <<ENDL \
#        set rd_use_mask \
#        set markersize 1 \
#        set grid vertical \
#        set output_filename {output.png} \
#        manhattan \
#        circular \
#        rd \
#        likelihood \
#        baf \
#        ENDL &&
#        cnvpytor -root {input.pytor} -view 100000 < script.spytor &&
#        &> {log}"""
