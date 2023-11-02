__author__ = "Jessika Nordin, Patrik Smeds"
__copyright__ = "Copyright 2022, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvpytor_readdepth:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        vcf="{sample}_{type}.germline.vcf",
    output:
        pytor=temp("cnv_sv/cnvpytor/{sample}_{type}.pytor"),
    params:
        extra=config.get("cnvpytor_readdepth", {}).get("extra", ""),
        length=config.get("cnvpytor_readdepth", {}).get("length_list", ""),
        model=config.get("cnvpytor_readdepth", {}).get("calling_model", ""),
    log:
        "cnv_sv/cnvpytor/{sample}_{type}_rd.log",
    benchmark:
        repeat("cnv_sv/cnvpytor/{sample}_{type}_rd.benchmark.tsv", config.get("cnvpytor_readdepth", {}).get("benchmark_repeats", 1))
    threads: config.get("cnvpytor_readdepth", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvpytor_readdepth", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvpytor_readdepth", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvpytor_readdepth", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvpytor_readdepth", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvpytor_readdepth", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvpytor_readdepth", {}).get("container", config["default_container"])
    message:
        "{rule}: Run cnvpytor calls based on read depth for {wildcards.sample}_{wildcards.type}"
    shell:
        """cnvpytor --max_cores {threads} -root {output.pytor} {params.extra} -rd {input.bam} &> {log} &&
        cnvpytor --max_cores {threads}  -root {output.pytor} -his {params.length} &>> {log} &&
        cnvpytor --max_cores {threads}  -root {output.pytor} -partition {params.length} &>> {log} &&
        cnvpytor --max_cores {threads}  -root {output.pytor} -snp {input.vcf} -sample {wildcards.sample}_{wildcards.type} \
         &>> {log} &&
        cnvpytor --max_cores {threads}  -root {output.pytor} -mask_snps &>> {log} &&
        cnvpytor --max_cores {threads}  -root {output.pytor} -baf {params.length} &>> {log} &&
        cnvpytor --max_cores {threads}  -root {output.pytor} -call {params.model} {params.length} &>> {log}
        """


rule cnvpytor_filter:
    input:
        pytor="cnv_sv/cnvpytor/{sample}_{type}.pytor",
    output:
        vcf=temp("cnv_sv/cnvpytor/{sample}_{type}.vcf"),
        filtvcf=temp("cnv_sv/cnvpytor/{sample}_{type}.filtered.vcf"),
    params:
        extra=config.get("cnvpytor_filter", {}).get("extra", ""),
        dgrange=config.get("cnvpytor_filter", {}).get("dG_range", ""),
        model=config.get("cnvpytor_filter", {}).get("calling_model", "rd_mean_shift"),
        prange=config.get("cnvpytor_filter", {}).get("p_range", ""),
        pnrange=config.get("cnvpytor_filter", {}).get("pN_range", ""),
        q0range=config.get("cnvpytor_filter", {}).get("Q0_range", ""),
        view=config.get("cnvpytor_filter", {}).get("view", ""),
    log:
        "cnv_sv/cnvpytor/{sample}_{type}_filter.log",
    benchmark:
        repeat("cnv_sv/cnvpytor/{sample}_{type}_filter.benchmark.tsv", config.get("cnvpytor_filter", {}).get("benchmark_repeats", 1))
    threads: config.get("cnvpytor_filter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvpytor_filter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvpytor_filter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvpytor_filter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvpytor_filter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvpytor_filter", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvpytor_filter", {}).get("container", config["default_container"])
    message:
        "{rule}: Filter cnvpytor calls for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        cnvpytor -root {input.pytor} -view {params.view} <<-ENDL  &> {log}
        set callers {params.model}
        set print_filename {output.vcf}
        print calls
        set Q0_range {params.q0range}
        set p_range {params.prange}
        set pN_range {params.pnrange}
        set dG_range {params.dgrange}
        set print_filename {output.filtvcf}
        print calls
        ENDL
        """


## Might be added later, code not tested ##
# rule cnvpytor_plot:
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
