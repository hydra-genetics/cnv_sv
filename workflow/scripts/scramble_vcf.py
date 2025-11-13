from datetime import date

meis_in = open(snakemake.input.meis)
vcf_out = open(snakemake.output.vcf, "w")
sample_name = snakemake.params.sample_name
caller = snakemake.params.caller


def write_vcf_header(sample_name):
    vcf_out.write("##fileformat=VCFv4.2\n")
    vcf_out.write("##fileDate=%s\n" % str(date.today()))
    vcf_out.write("##source=SCRAMBLE\n")
    vcf_out.write(
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
    vcf_out.write(
        "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    vcf_out.write(
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    vcf_out.write(
        "##INFO=<ID=MEINFO,Number=4,Type=String,Description=\"Mobile element info of the form NAME,START,END,POLARITY\">\n")
    vcf_out.write(
        "##INFO=<ID=CALLER,Number=1,Type=String,Description=\"Caller\">\n")
    vcf_out.write(
        "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of supporting reads\">\n")
    vcf_out.write(
        "##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"SCRAMBLE alignment score\">\n")
    vcf_out.write(
        "##INFO=<ID=CLIP_SIDE,Number=1,Type=String,Description=\"Side of soft-clipping (left/right)\">\n")
    vcf_out.write(
        "##INFO=<ID=POLYA_LEN,Number=1,Type=Integer,Description=\"Length of polyA tail\">\n")
    vcf_out.write(
        "##INFO=<ID=POLYA_SCORE,Number=1,Type=Float,Description=\"PolyA alignment score\">\n")
    vcf_out.write(
        "##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">\n")
    vcf_out.write(
        "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">\n")
    vcf_out.write(
        "##ALT=<ID=INS:ME:SVA,Description=\"Insertion of SVA element\">\n")
    vcf_out.write(
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    vcf_out.write(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample_name)


header_map = {}
for line in meis_in:
    columns = line.strip().split("\t")
    if columns[0] == "Insertion" or columns[0].startswith("#"):
        header_map = {column_name: index for index,
                      column_name in enumerate(columns)}
        write_vcf_header(sample_name)
        continue
    if not line.strip():
        continue
    location = columns[header_map.get('Insertion', 0)]
    if ':' not in location:
        continue
    chrom, pos = location.split(':')
    try:
        pos_int = int(pos)
        if pos_int <= 0:
            continue
    except ValueError:
        continue
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    mei_type = columns[header_map.get('MEI_Family', 1)].upper()
    if mei_type == "LINE1":
        mei_type = "L1"
    orientation = columns[header_map.get('Orientation', 2)]
    polarity = "+" if orientation == "Plus" else "-"
    support = columns[header_map.get('Support', 3)]
    score = columns[header_map.get('Score', 4)]
    polya_len = columns[header_map.get('polyA_Position', 5)]
    polya_score = columns[header_map.get('polyA_Percent_Match', 6)]
    consensus = columns[header_map.get('Consensus', 7)]
    clip_side = columns[header_map.get('Clipped_Side', 8)]
    ref = "N"
    alt = f"<INS:ME:{mei_type}>"
    id = "."
    qual = "."
    filter = "PASS"

    # Estimate SVLEN (TODO: add actual consensus length if available)
    svlen = "."
    if mei_type == "ALU":
        svlen = "300"
    elif mei_type == "L1":
        svlen = "6000"
    elif mei_type == "SVA":
        svlen = "2000"

    end = str(pos_int + 1)
    info = f"SVTYPE=INS;END={end}"
    if svlen != ".":
        info += f";SVLEN={svlen}"
    info += f";MEINFO={mei_type},{pos},{end},{polarity}"
    info += f";CALLER={caller}"
    info += f";SUPPORT={support}"
    if score and score != "NA":
        info += f";SCORE={score}"
    if clip_side and clip_side != "NA":
        info += f";CLIP_SIDE={clip_side}"
    if polya_len and polya_len != "NA" and polya_len != "None Found":
        info += f";POLYA_LEN={polya_len}"
    if polya_score and polya_score != "NA" and polya_score != "None Found":
        info += f";POLYA_SCORE={polya_score}"

    # Genotype - assume heterozygous for MEI insertions (TODO: doublecheck assumption)
    format_field = "GT"
    data = "0/1"

    out_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
        chrom, pos, id, ref, alt, qual, filter, info, format_field, data
    )
    vcf_out.write(out_line)
vcf_out.close()
