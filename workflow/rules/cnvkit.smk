__author__ = "Jonas Almlöf, Martin Rippin"
__copyright__ = "Copyright 2022, Jonas Almlöf, Martin Rippin"
__email__ = "jonas.almlof@scilifelab.uu.se, martin.rippin@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_batch:
    input:
        bam=lambda wildcards: get_haplotagged_input_bam(wildcards)[0],
        bai=lambda wildcards: get_haplotagged_input_bam(wildcards)[1],
        reference=config.get("cnvkit_batch", {}).get("normal_reference", ""),
    output:
        antitarget_coverage=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.antitargetcoverage.cnn"),
        bins=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.bintest.cns"),
        regions=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr"),
        segments=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns"),
        segments_called=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.call.cns"),
        target_coverage=temp("cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.targetcoverage.cnn"),
    params:
        extra=config.get("cnvkit_batch", {}).get("extra", ""),
        method=config.get("cnvkit_batch", {}).get("method", "hybrid"),
    log:
        "cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.call.cns.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.call.cns.benchmark.tsv",
            config.get("cnvkit_batch", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_batch", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_batch", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_batch", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_batch", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_batch", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_batch", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit_batch", {}).get("container", config["default_container"])
    message:
        "{rule}: use cnvkit to call cnvs in {wildcards.sample}/{wildcards.sample}_{wildcards.type}"
    wrapper:
        "v3.3.6/bio/cnvkit/batch"


rule cnvkit_call:
    input:
        segment="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.germline.vcf",
        tc_file=get_tc_file,
    output:
        segment=temp("cnv_sv/cnvkit_call/{sample}_{type}.{tc_method}.loh.cns"),
    params:
        extra=config.get("cnvkit_call", {}).get("extra", ""),
        purity=get_tc,
    log:
        "cnv_sv/cnvkit_call/{sample}_{type}.{tc_method}.loh.cns.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_call/{sample}_{type}.{tc_method}.loh.cns.benchmark.tsv",
            config.get("cnvkit_call", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_call", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_call", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit_call", {}).get("container", config["default_container"])
    message:
        "{rule}: call cnvs with loh info into cnv_sv/cnvkit_call/{wildcards.sample}_{wildcards.type}.loh.cns"
    wrapper:
        "v3.3.6/bio/cnvkit/call"


rule cnvkit_diagram:
    input:
        cns="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
        cnr="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr",
    output:
        pdf=temp("cnv_sv/cnvkit_diagram/{sample}_{type}.pdf"),
    params:
        extra=config.get("cnvkit_diagram", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_diagram/{sample}_{type}.pdf.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_diagram/{sample}_{type}.pdf.benchmark.tsv",
            config.get("cnvkit_diagram", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_diagram", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_diagram", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_diagram", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_diagram", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_diagram", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_diagram", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit_diagram", {}).get("container", config["default_container"])
    message:
        "{rule}: chromosome plot cnv_sv/cnvkit_scatter/{wildcards.sample}_{wildcards.type}.pdf"
    wrapper:
        "v3.3.6/bio/cnvkit/diagram"


rule cnvkit_scatter:
    input:
        segments="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
        segment_regions="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr",
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.germline.vcf",
    output:
        plot=temp("cnv_sv/cnvkit_scatter/{sample}_{type}.png"),
    params:
        extra=config.get("cnvkit_scatter", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_scatter/{sample}_{type}.png.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_scatter/{sample}_{type}.png.benchmark.tsv",
            config.get("cnvkit_scatter", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_scatter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_scatter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_scatter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_scatter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_scatter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_scatter", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit_scatter", {}).get("container", config["default_container"])
    message:
        "{rule}: plot cnvs into cnv_sv/cnvkit_scatter/{wildcards.sample}_{wildcards.type}.png"
    shell:
        "(cnvkit.py scatter {input.segment_regions} "
        "-s {input.segments} "
        "-v {input.vcf} "
        "-o {output.plot} "
        "{params.extra}) &> {log}"


rule cnvkit_vcf:
    input:
        segment="cnv_sv/cnvkit_call/{sample}_{type}.{tc_method}.loh.cns",
    output:
        vcf=temp("cnv_sv/cnvkit_vcf/{sample}_{type}.{tc_method}.vcf"),
    params:
        sample_name="{sample}_{type}",
        hom_del_limit=config.get("cnvkit_vcf", {}).get("hom_del_limit", 0.5),
        het_del_limit=config.get("cnvkit_vcf", {}).get("het_del_limit", 1.5),
        dup_limit=config.get("cnvkit_vcf", {}).get("dup_limit", 2.5),
        caller=config.get("cnvkit_vcf", {}).get("caller_name", "cnvkit"),
    log:
        "cnv_sv/cnvkit_vcf/{sample}_{type}.{tc_method}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_vcf/{sample}_{type}.{tc_method}.vcf.benchmark.tsv",
            config.get("cnvkit_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit_vcf", {}).get("container", config["default_container"])
    message:
        "{rule}: export cnvkit segments into vcf in {output.vcf}"
    script:
        "../scripts/cnvkit_vcf.py"


rule cnvkit_export_seg:
    input:
        segments="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
    output:
        seg=temp("cnv_sv/cnvkit_export_seg/{sample}_{type}.seg"),
    params:
        extra=config.get("cnvkit_export_seg", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_export_seg/{sample}_{type}.seg.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_export_seg/{sample}_{type}.seg.benchmark.tsv",
            config.get("cnvkit_export_seg", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_export_seg", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_export_seg", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_export_seg", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_export_seg", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_export_seg", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_export_seg", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit_export_seg", {}).get("container", config["default_container"])
    message:
        "{rule}: export cnvkit segments into seg in {output.seg}"
    wrapper:
        "v3.3.6/bio/cnvkit/export"
