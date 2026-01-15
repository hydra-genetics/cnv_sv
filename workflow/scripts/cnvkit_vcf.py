import logging
import math
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


def write_vcf_header(
    vcf_out: Any, sample_name: str, caller: str, date_str: Optional[str] = None
) -> None:
    """Writes the VCF header to the output file.

    Args:
        vcf_out: Opened file handle for VCF output.
        sample_name: Name of the sample.
        caller: Name of the CNV caller.
        date_str: Optional date string. Defaults to today's date.
    """
    if date_str is None:
        date_str = str(date.today())

    vcf_out.write("##fileformat=VCFv4.2\n")
    vcf_out.write(f"##fileDate={date_str}\n")
    vcf_out.write(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n'
    )
    vcf_out.write(
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'
    )
    vcf_out.write(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
    )
    vcf_out.write(f'##INFO=<ID=CALLER,Number=1,Type=String,Description="Caller={caller}">\n')
    vcf_out.write('##INFO=<ID=CN,Number=1,Type=Float,Description="Copy number">\n')
    vcf_out.write(
        '##INFO=<ID=CORR_CN,Number=1,Type=Float,Description="Tumour content corrected copy number">\n'
    )
    vcf_out.write(
        '##INFO=<ID=LOG_ODDS_RATIO,Number=1,Type=Float,Description="Log odds ratio">\n'
    )
    vcf_out.write(
        '##INFO=<ID=PROBES,Number=1,Type=Integer,Description="Number of probes in CNV">\n'
    )
    vcf_out.write('##INFO=<ID=BAF,Number=1,Type=Float,Description="SNP minor allele frequency">\n')
    vcf_out.write('##ALT=<ID=DEL,Description="Deletion">\n')
    vcf_out.write('##ALT=<ID=DUP,Description="Duplication">\n')
    vcf_out.write('##ALT=<ID=COPY_NORMAL,Description="Normal copy number">\n')
    vcf_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf_out.write('##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number">\n')
    vcf_out.write('##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Number of probes in CNV">\n')
    vcf_out.write(
        '##FORMAT=<ID=DP,Number=1,Type=Float,Description="Average coverage over region">\n'
    )
    vcf_out.write('##FORMAT=<ID=BAF,Number=1,Type=Float,Description="SNP minor allele frequency">\n')
    vcf_out.write(
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
    )


def calculate_cn_and_gt(
    log2: float, hom_del_limit: float, het_del_limit: float, dup_limit: float
) -> Tuple[float, str, str]:
    """Calculates copy number, SV type, and genotype from log2 ratio.

    Args:
        log2: Log2 ratio from CNVkit.
        hom_del_limit: Threshold for homozygous deletion.
        het_del_limit: Threshold for heterozygous deletion.
        dup_limit: Threshold for duplication.

    Returns:
        Tuple containing (copy number, SV type, genotype).
    """
    cn = round(2 * pow(2, log2), 2)

    if cn < het_del_limit:
        alt = "<DEL>"
    elif cn > dup_limit:
        alt = "<DUP>"
    else:
        alt = "<COPY_NORMAL>"

    if cn < hom_del_limit:
        gt = "1/1"
    elif (cn >= hom_del_limit and cn < het_del_limit) or cn > dup_limit:
        gt = "0/1"
    else:
        gt = "0/0"

    return cn, alt, gt


def create_vcf_record(
    chrom: str,
    start: str,
    end: str,
    log2: str,
    probes: str,
    baf: str,
    depth: str,
    caller: str,
    cn: float,
    alt: str,
    gt: str,
) -> str:
    """Creates a VCF record line.

    Args:
        chrom: Chromosome name.
        start: Start position.
        end: End position.
        log2: Log2 ratio.
        probes: Number of probes.
        baf: B-allele frequency.
        depth: Average depth.
        caller: Caller name.
        cn: Calculated copy number.
        alt: SV type (e.g., <DEL>, <DUP>).
        gt: Genotype.

    Returns:
        VCF record as a string.
    """
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"

    svlen = int(end) - int(start) + 1
    ref = "N"
    qual = "."
    vcf_filter = "."
    vcf_id = "."

    # INFO field
    info_items = [
        f"SVTYPE={alt[1:-1]}",
        f"END={end}",
        f"SVLEN={svlen}",
        f"LOG_ODDS_RATIO={log2}",
        f"CALLER={caller}",
        "CN=NA",
        f"CORR_CN={cn}",
        f"PROBES={probes}",
        f"BAF={baf}" if baf else f"BAF={0.0}",
    ]
    info = ";".join(info_items)

    # FORMAT and Sample data
    format_fields = ["GT", "CN", "CNQ", "DP"]
    data_fields = [gt, str(cn), probes, depth]

    if baf:
        format_fields.append("BAF")
        data_fields.append(baf)

    vcf_format = ":".join(format_fields)
    vcf_data = ":".join(data_fields)

    return f"{chrom}\t{start}\t{vcf_id}\t{ref}\t{alt}\t{qual}\t{vcf_filter}\t{info}\t{vcf_format}\t{vcf_data}\n"


def process_segments(
    seg_path: str,
    vcf_path: str,
    sample_name: str,
    caller: str,
    hom_del_limit: float,
    het_del_limit: float,
    dup_limit: float,
) -> None:
    """Processes CNVkit segment file and writes to VCF.

    Args:
        seg_path: Path to the input segment file.
        vcf_path: Path to the output VCF file.
        sample_name: Name of the sample.
        caller: Caller name.
        hom_del_limit: Threshold for homozygous deletion.
        het_del_limit: Threshold for heterozygous deletion.
        dup_limit: Threshold for duplication.
    """
    header_map: Dict[str, int] = {}
    with open(seg_path, "r") as seg_in, open(vcf_path, "w") as vcf_out:
        for line in seg_in:
            columns = line.strip().split("\t")
            if not columns:
                continue

            if columns[0] == "chromosome":
                header_map = {name: i for i, name in enumerate(columns)}
                write_vcf_header(vcf_out, sample_name, caller)
                continue

            if not header_map:
                logging.warning("Segment file missing header. Skipping line: %s", line)
                continue

            try:
                chrom = columns[header_map["chromosome"]]
                start = columns[header_map["start"]]
                end = columns[header_map["end"]]
                probes = columns[header_map["probes"]]
                log2 = columns[header_map["log2"]]
                baf = columns[header_map.get("baf", "")] if "baf" in header_map else ""
                depth = columns[header_map.get("depth", "")] if "depth" in header_map else "0"

                cn, alt, gt = calculate_cn_and_gt(
                    float(log2), hom_del_limit, het_del_limit, dup_limit
                )

                vcf_record = create_vcf_record(
                    chrom, start, end, log2, probes, baf, depth, caller, cn, alt, gt
                )
                vcf_out.write(vcf_record)
            except (KeyError, IndexError, ValueError) as e:
                logging.error("Error processing line: %s. Error: %s", line.strip(), e)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    process_segments(
        seg_path=snakemake.input.segment,
        vcf_path=snakemake.output.vcf,
        sample_name=snakemake.params.sample_name,
        caller=snakemake.params.caller,
        hom_del_limit=snakemake.params.hom_del_limit,
        het_del_limit=snakemake.params.het_del_limit,
        dup_limit=snakemake.params.dup_limit,
    )
