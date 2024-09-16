
from datetime import date

seg_in = open(snakemake.input.segment)
vcf_out = open(snakemake.output.vcf, "w")
sample_name = snakemake.params.sample_name
hom_del_limit = snakemake.params.hom_del_limit
het_del_limit = snakemake.params.het_del_limit
dup_limit = snakemake.params.dup_limit


def write_vcf_header(gatk_version, sample_name):
    vcf_out.write("##fileformat=VCFv4.2\n")
    vcf_out.write("##fileDate=%s\n" % str(date.today()))
    vcf_out.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
    vcf_out.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    vcf_out.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    vcf_out.write("##INFO=<ID=CALLER,Number=1,Type=String,Description=\"Caller\">\n")
    vcf_out.write("##INFO=<ID=CN,Number=1,Type=Float,Description=\"Copy number\">\n")
    vcf_out.write("##INFO=<ID=CORR_CN,Number=1,Type=Float,Description=\"Tumour content corrected copy number\">\n")
    vcf_out.write("##INFO=<ID=LOG_ODDS_RATIO,Number=1,Type=Float,Description=\"Log odds ratio\">\n")
    vcf_out.write("##INFO=<ID=PROBES,Number=1,Type=Integer,Description=\"Number of probes in CNV\">\n")
    vcf_out.write("##INFO=<ID=BAF,Number=1,Type=Float,Description=\"SNP minor allele frequency\">\n")
    vcf_out.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
    vcf_out.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
    vcf_out.write("##ALT=<ID=COPY_NORMAL,Description=\"Normal copy number\">\n")
    vcf_out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    vcf_out.write("##FORMAT=<ID=CN,Number=1,Type=Float,Description=\"Copy number\">\n")
    vcf_out.write("##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"Number of probes in CNV\">\n")
    vcf_out.write("##FORMAT=<ID=DP,Number=1,Type=Float,Description=\"Average coverage over region\">\n")
    vcf_out.write("##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"SNP minor allele frequency\">\n")
    vcf_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample_name)


header = 0
gatk_version = ""
header_map = {}
for line in seg_in:
    columns = line.strip().split("\t")
    if columns[0] == "chromosome":
        header_map = {column_name: index for index, column_name in enumerate(columns)}
        write_vcf_header(gatk_version, sample_name)
    else:
        chrom = columns[header_map['chromosome']]
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        start_pos = columns[header_map['start']]
        end_pos = columns[header_map['end']]
        svlen = int(end_pos) - int(start_pos) + 1
        nr_probes = columns[header_map['probes']]
        log_odds_ratio = columns[header_map['log2']]
        baf = columns[header_map['baf']]
        dp = columns[header_map['depth']]
        cn = round(2*pow(2, float(log_odds_ratio)), 2)
        ref = "N"
        alt = ""
        if cn < het_del_limit:
            alt = "<DEL>"
        elif cn > dup_limit:
            alt = "<DUP>"
        else:
            alt = "<COPY_NORMAL>"
        id = "."
        qual = "."
        filter = "."
        gt = ""
        if cn < hom_del_limit:
            gt = "1/1"
        elif (cn >= hom_del_limit and cn < het_del_limit) or cn > dup_limit:
            gt = "0/1"
        else:
            gt = "0/0"
        info = "SVTYPE=%s;END=%s;SVLEN=%s;LOG_ODDS_RATIO=%s;CALLER=cnvkit;CN=NA;CORR_CN=%s" % (
            alt[1:-1], end_pos, svlen, log_odds_ratio, str(cn)
        )
        info = "%s;PROBES=%s;BAF=%s" % (info, nr_probes, baf)
        format = "GT:CN:CNQ:DP"
        data = "%s:%s:%s:%s" % (gt, cn, nr_probes, dp)
        if baf != "":
            format += ":BAF"
            data += ":%s" % (baf)
        out_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, start_pos, id, ref, alt, qual, filter, info, format, data)
        vcf_out.write(out_line)
vcf_out.close()
