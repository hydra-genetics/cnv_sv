
from datetime import date

seg_in = open(snakemake.input.segment)
vcf_out = open(snakemake.output.vcf, "w")
hom_del_limit = snakemake.params.hom_del_limit
het_del_limit = snakemake.params.het_del_limit
dup_limit = snakemake.params.dup_limit


def write_vcf_header(sample_name):
    vcf_out.write("##fileformat=VCFv4.2\n")
    vcf_out.write("##fileDate=%s\n" % str(date.today()))
    vcf_out.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
    vcf_out.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    vcf_out.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    vcf_out.write("##INFO=<ID=COPY_NUMBER,Number=1,Type=Float,Description=\"Copy number\">\n")
    vcf_out.write("##INFO=<ID=LOG_ODDS_RATIO,Number=1,Type=Float,Description=\"Log odds ratio\">\n")
    vcf_out.write("##INFO=<ID=PROBES,Number=1,Type=Integer,Description=\"Number of probes in CNV\">\n")
    vcf_out.write("##INFO=<ID=BAF,Number=1,Type=Float,Description=\"SNP minor allele frequency\">\n")
    vcf_out.write("##INFO=<ID=BAF_PROBES,Number=1,Type=Integer,Description=\"Number of SNPs for BAF\">\n")
    vcf_out.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
    vcf_out.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
    vcf_out.write("##ALT=<ID=COPY_NORMAL,Description=\"Normal copy number\">\n")
    vcf_out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    vcf_out.write("##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">\n")
    vcf_out.write("##FORMAT=<ID=CNQ,Number=1,Type=Integer,Description=\"Number of probes in CNV\">\n")
    vcf_out.write("##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"SNP minor allele frequency\">\n")
    vcf_out.write("##FORMAT=<ID=BAFQ,Number=1,Type=Integer,Description=\"Number of SNPs for BAF\">\n")
    vcf_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample_name)


header = 0
gatk_version = ""
sample_name = ""
for line in seg_in:
    columns = line.strip().split("\t")
    if columns[0] == "@RG":
        sample_name = columns[2].split(":")[1]
    elif columns[0] == "CONTIG":
        write_vcf_header(sample_name)
    elif columns[0][0] == "@":
        continue
    else:
        chrom = columns[0]
        start_pos = columns[1]
        end_pos = columns[2]
        svlen = int(end_pos) - int(start_pos) + 1
        nr_probes = columns[3]
        log_odds_ratio = columns[6]
        nr_baf_probes = columns[4]
        baf = columns[9]
        cn = round(2*pow(2, float(log_odds_ratio)), 1)
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
        info = "SVTYPE=%s;END=%s;SVLEN=%s;LOG_ODDS_RATIO=%s;COPY_NUMBER=%s;PROBES=%s;BAF=%s;BAF_PROBES=%s" % (
            ref, end_pos, svlen, log_odds_ratio, str(cn), nr_probes, baf, nr_baf_probes
        )
        format = "GT:CN:CNQ"
        data = "%s:%s:%s" % (gt, cn, nr_probes)
        if int(nr_baf_probes) > 0:
            format += ":BAF:BAFQ"
            data += ":%s:%s" % (baf, nr_baf_probes)
        out_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, start_pos, id, ref, alt, qual, filter, info, format, data)
        vcf_out.write(out_line)
vcf_out.close()
