__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2026, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule ichorcna_offtarget_read_counter:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
    output:
        wig=temp("cnv_sv/ichorcna_offtarget_read_counter/{sample}_{type}.wig"),
    params:
        chrs=config.get("ichorcna_offtarget_read_counter", {}).get(
            "chrs",
            "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,"
            "chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX",
        ),
        window=config.get("ichorcna_offtarget_read_counter", {}).get("window", 100000),
        quality=config.get("ichorcna_offtarget_read_counter", {}).get("quality", 20),
    log:
        "cnv_sv/ichorcna_offtarget_read_counter/{sample}_{type}.wig.log",
    benchmark:
        repeat(
            "cnv_sv/ichorcna_offtarget_read_counter/{sample}_{type}.wig.benchmark.tsv",
            config.get("ichorcna_offtarget_read_counter", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("ichorcna_offtarget_read_counter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("ichorcna_offtarget_read_counter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("ichorcna_offtarget_read_counter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("ichorcna_offtarget_read_counter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("ichorcna_offtarget_read_counter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("ichorcna_offtarget_read_counter", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("ichorcna_offtarget_read_counter", {}).get("container", config["default_container"])
    message:
        "{rule}: Count reads in {params.window}bp bins for {input.bam} with HMMcopy readCounter"
    shell:
        "readCounter {input.bam} -c {params.chrs} -w {params.window} -q {params.quality} > {output.wig} 2> {log}"


rule ichorcna_offtarget_run:
    input:
        wig="cnv_sv/ichorcna_offtarget_read_counter/{sample}_{type}.wig",
    output:
        cna=temp("cnv_sv/ichorcna_offtarget_run/{sample}_{type}/{sample}_{type}.cna.seg"),
        seg=temp("cnv_sv/ichorcna_offtarget_run/{sample}_{type}/{sample}_{type}.seg"),
        seg_txt=temp("cnv_sv/ichorcna_offtarget_run/{sample}_{type}/{sample}_{type}.seg.txt"),
        params_txt=temp("cnv_sv/ichorcna_offtarget_run/{sample}_{type}/{sample}_{type}.params.txt"),
        corrected_depth=temp("cnv_sv/ichorcna_offtarget_run/{sample}_{type}/{sample}_{type}.correctedDepth.txt"),
        rdata=temp("cnv_sv/ichorcna_offtarget_run/{sample}_{type}/{sample}_{type}.RData"),
        plot_dir=temp(directory("cnv_sv/ichorcna_offtarget_run/{sample}_{type}/{sample}_{type}/")),
    params:
        out_dir="cnv_sv/ichorcna_offtarget_run/{sample}_{type}/",
        id="{sample}_{type}",
        # gc_wig/map_wig/centromere are params, not input: they're often paths
        # bundled inside the container (e.g. /opt/ichorCNA/inst/extdata/...),
        # which don't exist on the host filesystem Snakemake itself checks
        # against - declaring them as input would make dry-run/DAG-building
        # fail with a false "missing input file" for any such path.
        gc_wig=config.get("ichorcna_offtarget_run", {}).get("gc_wig", ""),
        map_wig=config.get("ichorcna_offtarget_run", {}).get("map_wig", ""),
        centromere=config.get("ichorcna_offtarget_run", {}).get("centromere", ""),
        ploidy=config.get("ichorcna_offtarget_run", {}).get("ploidy", "c(2,3,4)"),
        normal=config.get("ichorcna_offtarget_run", {}).get("normal", "c(0.5)"),
        max_cn=config.get("ichorcna_offtarget_run", {}).get("max_cn", 7),
        include_homd=config.get("ichorcna_offtarget_run", {}).get("include_homd", "FALSE"),
        chrs=config.get("ichorcna_offtarget_run", {}).get("chrs", 'c(1:22,"X")'),
        chr_train=config.get("ichorcna_offtarget_run", {}).get("chr_train", "c(1:22)"),
        genome_build=config.get("ichorcna_offtarget_run", {}).get("genome_build", "hg38"),
        genome_style=config.get("ichorcna_offtarget_run", {}).get("genome_style", "UCSC"),
        min_map_score=config.get("ichorcna_offtarget_run", {}).get("min_map_score", 0.9),
        extra=get_ichorcna_offtarget_run_extra,
    log:
        "cnv_sv/ichorcna_offtarget_run/{sample}_{type}/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/ichorcna_offtarget_run/{sample}_{type}/{sample}_{type}.output.benchmark.tsv",
            config.get("ichorcna_offtarget_run", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("ichorcna_offtarget_run", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("ichorcna_offtarget_run", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("ichorcna_offtarget_run", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("ichorcna_offtarget_run", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("ichorcna_offtarget_run", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("ichorcna_offtarget_run", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("ichorcna_offtarget_run", {}).get("container", config["default_container"])
    message:
        "{rule}: Estimate tumor fraction from off-target bins in {input.wig} with ichorCNA"
    shell:
        "(Rscript /opt/ichorCNA/scripts/runIchorCNA.R "
        "--id {params.id} "
        "--WIG {input.wig} "
        "--gcWig {params.gc_wig} "
        "--mapWig {params.map_wig} "
        "--centromere {params.centromere} "
        "--ploidy \"{params.ploidy}\" "
        "--normal \"{params.normal}\" "
        "--maxCN {params.max_cn} "
        "--includeHOMD {params.include_homd} "
        "--chrs '{params.chrs}' "
        "--chrTrain '{params.chr_train}' "
        "--genomeBuild {params.genome_build} "
        "--genomeStyle {params.genome_style} "
        "--minMapScore {params.min_map_score} "
        "{params.extra} "
        "--outDir {params.out_dir}) &> {log}"
