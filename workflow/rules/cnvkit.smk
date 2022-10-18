__author__ = "Jonas Almlöf, Martin Rippin"
__copyright__ = "Copyright 2022, Jonas Almlöf, Martin Rippin"
__email__ = "jonas.almlof@scilifelab.uu.se, martin.rippin@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_batch:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        cnv_reference=config.get("cnvkit_batch", {}).get("normal_reference", ""),
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
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
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
    conda:
        "../envs/cnvkit.yaml"
    message:
        "{rule}: use cnvkit to call cnvs in {wildcards.sample}/{wildcards.sample}_{wildcards.type}"
    shell:
        "(cnvkit.py batch {input.bam} "
        "-r {input.cnv_reference} "
        "-d {params.outdir} "
        "-m {params.method} "
        "{params.extra}) &> {log}"


rule cnvkit_call:
    input:
        segment="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.germline.vcf",
    output:
        segment=temp("cnv_sv/cnvkit_call/{sample}_{type}.loh.cns"),
    params:
        extra=config.get("cnvkit_call", {}).get("extra", ""),
        tc=lambda wildcards: get_sample(samples, wildcards)["tumor_content"],
    log:
        "cnv_sv/cnvkit_call/{sample}_{type}.loh.cns.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_call/{sample}_{type}.loh.cns.benchmark.tsv",
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
    conda:
        "../envs/cnvkit.yaml"
    message:
        "{rule}: call cnvs with loh info into cnv_sv/cnvkit_call/{wildcards.sample}_{wildcards.type}.loh.cns"
    shell:
        "(cnvkit.py call {input.segment} "
        "-v {input.vcf} "
        "-o {output.segment} "
        "--purity {params.tc} "
        "{params.extra}) &> {log}"


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
    conda:
        "../envs/cnvkit.yaml"
    message:
        "{rule}: chromosome plot cnv_sv/cnvkit_scatter/{wildcards.sample}_{wildcards.type}.pdf"
    shell:
        "(cnvkit.py diagram {input.cnr} "
        "-s {input.cns} "
        "-o {output.pdf} "
        "{params.extra}) &> {log}"


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
    conda:
        "../envs/cnvkit.yaml"
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
        segment="cnv_sv/cnvkit_call/{sample}_{type}.loh.cns",
    output:
        vcf=temp("cnv_sv/cnvkit_vcf/{sample}_{type}.vcf"),
    params:
        sample_name="{sample}_{type}",
        hom_del_limit=config.get("cnvkit_vcf", {}).get("hom_del_limit", 0.5),
        het_del_limit=config.get("cnvkit_vcf", {}).get("het_del_limit", 1.5),
        dup_limit=config.get("cnvkit_vcf", {}).get("dup_limit", 2.5),
    log:
        "cnv_sv/cnvkit_vcf/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_vcf/{sample}_{type}.vcf.benchmark.tsv",
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
    conda:
        "../envs/cnvkit.yaml"
    message:
        "{rule}: export cnvkit segments into vcf in cnv_sv/cnvkit_vcf/{wildcards.sample}_{wildcards.type}.vcf"
    script:
        "../scripts/cnvkit_vcf.py"


rule cnvkit_seg:
    input:
        segments="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
    output:
        seg=temp("cnv_sv/cnvkit_seg/{sample}_{type}.seg"),
    params:
        extra=config.get("cnvkit_seg", {}).get("extra", ""),
    log:
        "cnv_sv/cnvkit_seg/{sample}_{type}.seg.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_seg/{sample}_{type}.seg.benchmark.tsv",
            config.get("cnvkit_seg", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_seg", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_seg", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_seg", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_seg", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_seg", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_seg", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit_seg", {}).get("container", config["default_container"])
    conda:
        "../envs/cnvkit.yaml"
    message:
        "{rule}: export cnvkit segments into seg in cnv_sv/{rule}/{wildcards.sample}_{wildcards.type}.seg"
    shell:
        "(cnvkit.py export seg "
        "{input.segments} "
        "--enumerate-chroms "
        "-o {output.seg} "
        "{params.extra}) &> {log}"
