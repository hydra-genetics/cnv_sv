"""
Adapted from GATK-SV make_scramble_vcf.py
(https://github.com/broadinstitute/gatk-sv/blob/v1.1/src/sv-pipeline/scripts/make_scramble_vcf.py)

Creates a VCF from the Scramble output table. Clusters redundant calls and
estimates variant length based on MEI type and alignment information.
Applies quality filters (Alignment_Score, Clipped_Reads_In_Cluster,
Alignment_Percent_Length, Alignment_Percent_Identity).
"""

import argparse
from collections import deque
from datetime import date
import logging
import sys

import pandas as pd

log = logging.getLogger()

DEFAULT_MIN_ALIGNMENT_SCORE = 70.0
DEFAULT_MIN_CLIPPED_READS = 5
DEFAULT_MIN_ALIGNMENT_PERCENT_LENGTH = 90.0
DEFAULT_MIN_ALIGNMENT_PERCENT_IDENTITY = 90.0

DEFAULT_ALU_SIZE = 282
DEFAULT_SVA_SIZE = 1362
DEFAULT_L1_SIZE = 6023
DEFAULT_CLUSTER_DISTANCE = 300

FLOAT_COLUMNS = ["Alignment_Score", "Alignment_Percent_Length", "Alignment_Percent_Identity"]
INT_COLUMNS = [
    "Clipped_Reads_In_Cluster", "Start_In_MEI", "Stop_In_MEI",
    "polyA_Position", "polyA_SupportingReads", "TSD_length",
]

# When merging two-sided calls, keep all values for these columns
# (e.g. polyA_Position=123,456) since they may differ between the two
# break-ends. They are not numeric values that can be averaged.
MERGED_INT_COLUMNS = [
    "Start_In_MEI", "Stop_In_MEI",
    "polyA_Position", "polyA_SupportingReads", "TSD_length",
]


def preprocess_scramble_table(path):
    """Parse SCRAMble MEI table, splitting Insertion into chrom/pos columns.

    Args:
        path: Path to the SCRAMble MEI TSV file

    Returns:
        pandas DataFrame with chrom, pos, end columns added, or empty DataFrame
    """
    df = pd.read_table(path, sep="\t")
    if df.empty:
        return df

    df = df[df["Insertion"].str.contains(":", na=False)].copy()
    if df.empty:
        return df

    df[["chrom", "pos"]] = df["Insertion"].str.split(":", n=1, expand=True)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["pos"])
    df["pos"] = df["pos"].astype(int)
    df["end"] = df["pos"] + 1

    for col in FLOAT_COLUMNS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in INT_COLUMNS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def calculate_svlen(df, alu_size=DEFAULT_ALU_SIZE, sva_size=DEFAULT_SVA_SIZE,
                    l1_size=DEFAULT_L1_SIZE, cluster_distance=DEFAULT_CLUSTER_DISTANCE):
    """Calculate SVLEN using GATK-SV logic with call clustering.

    Performs clustering of nearby calls with opposing clipped sides and matching
    insertion direction. Estimates variant length from MEI alignment positions.
    See: https://github.com/broadinstitute/gatk-sv/discussions/734

    Args:
        df: DataFrame from preprocess_scramble_table
        alu_size: Size in bp of ALU elements
        sva_size: Size in bp of SVA elements
        l1_size: Size in bp of LINE1 elements
        cluster_distance: Maximum distance to collapse pairs of break-end calls

    Returns:
        DataFrame with svlen column added
    """
    if df.empty:
        return df

    family_lengths = {"l1": l1_size, "sva": sva_size, "alu": alu_size}

    def _svlen_one_sided(r):
        fam_len = family_lengths.get(r["MEI_Family"])
        if fam_len is None or pd.isna(r.get("Start_In_MEI")) or pd.isna(r.get("Stop_In_MEI")):
            return family_lengths.get(r["MEI_Family"], 0)
        if r["Clipped_Side"] == "right":
            if r["Insertion_Direction"] == "Minus":
                return int(r["Stop_In_MEI"])
            else:
                return fam_len - int(r["Start_In_MEI"])
        else:
            if r["Insertion_Direction"] == "Minus":
                return fam_len - int(r["Start_In_MEI"])
            else:
                return int(r["Stop_In_MEI"])

    def _svlen_two_sided(r1, r2):
        if r1["Clipped_Side"] == "right":
            r_right, r_left = r1, r2
        else:
            r_right, r_left = r2, r1
        if r_right["Insertion_Direction"] == "Minus":
            return int(r_right["Stop_In_MEI"]) - int(r_left["Start_In_MEI"])
        else:
            return int(r_left["Stop_In_MEI"]) - int(r_right["Start_In_MEI"])

    df = df.sort_values(["chrom", "pos"]).reset_index(drop=True)

    buffer = deque()
    results = []

    for _, row in df.iterrows():
        row_dict = row.to_dict()
        row_dict["sides"] = 1

        new_buffer = deque()
        found_match = False

        for item in buffer:
            if found_match:
                new_buffer.append(item)
            elif (item["chrom"] != row_dict["chrom"]
                  or abs(int(item["pos"]) - int(row_dict["pos"])) > cluster_distance):
                item["svlen"] = _svlen_one_sided(item)
                results.append(item)
            elif (item["Clipped_Side"] != row_dict["Clipped_Side"]
                  and item["Insertion_Direction"] == row_dict["Insertion_Direction"]):
                row_dict["pos"] = item["pos"]
                row_dict["svlen"] = _svlen_two_sided(item, row_dict)
                row_dict["sides"] = 2
                row_dict["Clipped_Reads_In_Cluster"] = (
                    int(item["Clipped_Reads_In_Cluster"]) + int(row_dict["Clipped_Reads_In_Cluster"])
                )
                for col in ("Alignment_Score", "Alignment_Percent_Length", "Alignment_Percent_Identity"):
                    row_dict[col] = (item[col] + row_dict[col]) / 2
                for col in MERGED_INT_COLUMNS:
                    vals = []
                    for r in (item, row_dict):
                        v = r.get(col)
                        if pd.notna(v):
                            vals.append(str(int(v)))
                    row_dict[col] = ",".join(vals) if vals else float("nan")
                log.info(
                    "Clustered %s:%s and %s:%s (distance=%d, MEI_Family=%s)",
                    item["chrom"], item["pos"],
                    row_dict["chrom"], row_dict["Insertion"].split(":")[1],
                    abs(int(item["pos"]) - int(row_dict["Insertion"].split(":")[1])),
                    row_dict["MEI_Family"],
                )
                results.append(row_dict)
                found_match = True
            else:
                new_buffer.appendleft(item)

        if not found_match:
            new_buffer.appendleft(row_dict)
        buffer = new_buffer

    for item in buffer:
        item["svlen"] = _svlen_one_sided(item)
        results.append(item)

    return pd.DataFrame(results)


def apply_filters(df, min_score, min_reads, min_pct_length, min_pct_identity):
    """Apply quality filters, setting FILTER to PASS or LowQual.

    Args:
        df: DataFrame with SCRAMble calls
        min_score: Minimum Alignment_Score for PASS
        min_reads: Minimum Clipped_Reads_In_Cluster for PASS
        min_pct_length: Minimum Alignment_Percent_Length for PASS
        min_pct_identity: Minimum Alignment_Percent_Identity for PASS

    Returns:
        DataFrame with FILTER column added
    """
    if df.empty:
        return df

    pass_mask = (
        (df["Alignment_Score"] >= min_score)
        & (df["Clipped_Reads_In_Cluster"] >= min_reads)
        & (df["Alignment_Percent_Length"] >= min_pct_length)
        & (df["Alignment_Percent_Identity"] >= min_pct_identity)
    )
    df["FILTER"] = "LowQual"
    df.loc[pass_mask, "FILTER"] = "PASS"
    for _, row in df[~pass_mask].iterrows():
        log.info(
            "LowQual %s:%s (score=%.1f, reads=%d, pct_len=%.1f, pct_id=%.1f)",
            row["chrom"], int(row["pos"]),
            row["Alignment_Score"], int(row["Clipped_Reads_In_Cluster"]),
            row["Alignment_Percent_Length"], row["Alignment_Percent_Identity"],
        )
    return df


def write_vcf_header(vcf_out, sample_name):
    """Write VCF header lines.

    Args:
        vcf_out: Output file handle
        sample_name: Sample name for the header line
    """
    vcf_out.write("##fileformat=VCFv4.2\n")
    vcf_out.write(f"##fileDate={date.today()}\n")
    vcf_out.write("##source=scramble\n")
    vcf_out.write('##ALT=<ID=INS,Description="Insertion">\n')
    vcf_out.write(
        '##INFO=<ID=END,Number=1,Type=Integer,'
        'Description="End position of the variant described in this record">\n')
    vcf_out.write(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,'
        'Description="Type of structural variant">\n')
    vcf_out.write(
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,'
        'Description="Difference in length between REF and ALT alleles">\n')
    vcf_out.write(
        '##INFO=<ID=MEINFO,Number=4,Type=String,'
        'Description="Mobile element info of the form NAME,START,END,POLARITY">\n')
    vcf_out.write(
        '##INFO=<ID=ALGORITHMS,Number=.,Type=String,'
        'Description="Source algorithms">\n')
    vcf_out.write(
        '##INFO=<ID=Alignment_Score,Number=1,Type=Float,'
        'Description="SCRAMble alignment score">\n')
    vcf_out.write(
        '##INFO=<ID=Alignment_Percent_Length,Number=1,Type=Float,'
        'Description="SCRAMble alignment percent length">\n')
    vcf_out.write(
        '##INFO=<ID=Alignment_Percent_Identity,Number=1,Type=Float,'
        'Description="SCRAMble alignment percent identity">\n')
    vcf_out.write(
        '##INFO=<ID=Clipped_Reads_In_Cluster,Number=1,Type=Integer,'
        'Description="Number of clipped reads in cluster">\n')
    vcf_out.write(
        '##INFO=<ID=Start_In_MEI,Number=.,Type=Integer,'
        'Description="Start of alignment to canonical MEI sequence">\n')
    vcf_out.write(
        '##INFO=<ID=Stop_In_MEI,Number=.,Type=Integer,'
        'Description="End of alignment to canonical MEI sequence">\n')
    vcf_out.write(
        '##INFO=<ID=polyA_Position,Number=.,Type=Integer,'
        'Description="Position of poly-A tail">\n')
    vcf_out.write(
        '##INFO=<ID=polyA_SupportingReads,Number=.,Type=Integer,'
        'Description="Number of reads supporting poly-A tail">\n')
    vcf_out.write(
        '##INFO=<ID=TSD_length,Number=.,Type=Integer,'
        'Description="Target site duplication length">\n')
    vcf_out.write(
        '##FILTER=<ID=LowQual,Description="Low quality variant">\n')
    vcf_out.write(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf_out.write(
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")


def write_vcf_records(vcf_out, df, sample_name):
    """Write VCF data records from processed DataFrame.

    Args:
        vcf_out: Output file handle
        df: Processed DataFrame with svlen and FILTER columns
        sample_name: Sample name for record IDs
    """
    if df.empty:
        return

    for idx, (_, row) in enumerate(df.iterrows()):
        chrom = row["chrom"]
        pos = int(row["pos"])
        record_id = f"scramble_{sample_name}_{chrom}_{idx}"
        filt = row.get("FILTER", ".")
        polarity = "+" if row.get("Insertion_Direction") == "Plus" else "-"
        svlen = int(row["svlen"]) if pd.notna(row.get("svlen")) else 0

        info_parts = [
            f"SVTYPE=INS",
            f"END={pos}",
            f"SVLEN={svlen}",
            f"MEINFO={row['MEI_Family']},{pos},{pos},{polarity}",
            "ALGORITHMS=scramble",
        ]
        for col in FLOAT_COLUMNS:
            if pd.notna(row.get(col)):
                info_parts.append(f"{col}={row[col]:.2f}")
        if pd.notna(row.get("Clipped_Reads_In_Cluster")):
            info_parts.append(f"Clipped_Reads_In_Cluster={int(row['Clipped_Reads_In_Cluster'])}")
        for col in MERGED_INT_COLUMNS:
            val = row.get(col)
            if isinstance(val, str):
                info_parts.append(f"{col}={val}")
            elif pd.notna(val):
                info_parts.append(f"{col}={int(val)}")

        info = ";".join(info_parts)
        vcf_out.write(
            f"{chrom}\t{pos}\t{record_id}\tN\t<INS>\t.\t{filt}\t{info}\tGT\t0/1\n")


def process_meis_to_vcf(
        meis_file, vcf_out, sample_name,
        min_score=DEFAULT_MIN_ALIGNMENT_SCORE,
        min_reads=DEFAULT_MIN_CLIPPED_READS,
        min_pct_length=DEFAULT_MIN_ALIGNMENT_PERCENT_LENGTH,
        min_pct_identity=DEFAULT_MIN_ALIGNMENT_PERCENT_IDENTITY,
        alu_size=DEFAULT_ALU_SIZE,
        sva_size=DEFAULT_SVA_SIZE,
        l1_size=DEFAULT_L1_SIZE,
        cluster_distance=DEFAULT_CLUSTER_DISTANCE):
    """Process SCRAMble MEI file and write to VCF format.

    Adapted from GATK-SV make_scramble_vcf.py. Reads the MEI table with pandas,
    clusters nearby calls, estimates SVLEN, applies quality filters, and writes VCF.

    Args:
        meis_file: Path to SCRAMble MEI TSV file
        vcf_out: Output file handle for VCF
        sample_name: Sample name
        min_score: Minimum Alignment_Score for PASS filter
        min_reads: Minimum Clipped_Reads_In_Cluster for PASS filter
        min_pct_length: Minimum Alignment_Percent_Length for PASS filter
        min_pct_identity: Minimum Alignment_Percent_Identity for PASS filter
        alu_size: Size in bp of ALU elements
        sva_size: Size in bp of SVA elements
        l1_size: Size in bp of LINE1 elements
        cluster_distance: Maximum distance to collapse pairs of break-end calls
    """
    write_vcf_header(vcf_out, sample_name)

    df = preprocess_scramble_table(meis_file)
    if df.empty:
        return

    df = calculate_svlen(df, alu_size, sva_size, l1_size, cluster_distance)
    df.to_csv("debug_svlens.csv", index=False)  # Debug output to check SVLEN calculations
    if df.empty:
        return

    df = apply_filters(df, min_score, min_reads, min_pct_length, min_pct_identity)
    write_vcf_records(vcf_out, df, sample_name)


def main():
    try:
        snakemake  # noqa: F821
        snakemake.log_fmt_shell(stdout=False, stderr=True)  # noqa: F821
        meis_file = snakemake.input.meis  # noqa: F821
        vcf_path = snakemake.output.vcf  # noqa: F821
        sample_name = snakemake.params.sample_name  # noqa: F821
        min_score = float(snakemake.params.min_alignment_score)  # noqa: F821
        min_reads = int(snakemake.params.min_clipped_reads)  # noqa: F821
        min_pct_length = float(snakemake.params.min_alignment_percent_length)  # noqa: F821
        min_pct_identity = float(snakemake.params.min_alignment_percent_identity)  # noqa: F821
        alu_size = int(snakemake.params.get("alu_size", DEFAULT_ALU_SIZE))  # noqa: F821
        sva_size = int(snakemake.params.get("sva_size", DEFAULT_SVA_SIZE))  # noqa: F821
        l1_size = int(snakemake.params.get("l1_size", DEFAULT_L1_SIZE))  # noqa: F821
        cluster_distance = int(snakemake.params.get("cluster_distance", DEFAULT_CLUSTER_DISTANCE))  # noqa: F821
    except NameError:
        parser = argparse.ArgumentParser(
            description="Convert SCRAMble MEI output to VCF format. "
                        "Adapted from GATK-SV make_scramble_vcf.py.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        parser.add_argument("--input", required=True,
                            help="Input SCRAMble MEI file")
        parser.add_argument("--output", required=True,
                            help="Output VCF file")
        parser.add_argument("--sample-name", required=True,
                            help="Sample name")
        parser.add_argument("--min-alignment-score", type=float,
                            default=DEFAULT_MIN_ALIGNMENT_SCORE,
                            help="Minimum alignment score for PASS filter")
        parser.add_argument("--min-clipped-reads", type=int,
                            default=DEFAULT_MIN_CLIPPED_READS,
                            help="Minimum clipped reads for PASS filter")
        parser.add_argument("--min-alignment-percent-length", type=float,
                            default=DEFAULT_MIN_ALIGNMENT_PERCENT_LENGTH,
                            help="Minimum alignment percent length for PASS filter")
        parser.add_argument("--min-alignment-percent-identity", type=float,
                            default=DEFAULT_MIN_ALIGNMENT_PERCENT_IDENTITY,
                            help="Minimum alignment percent identity for PASS filter")
        parser.add_argument("--alu-size", type=int, default=DEFAULT_ALU_SIZE,
                            help="ALU element size in bp")
        parser.add_argument("--sva-size", type=int, default=DEFAULT_SVA_SIZE,
                            help="SVA element size in bp")
        parser.add_argument("--l1-size", type=int, default=DEFAULT_L1_SIZE,
                            help="LINE1 element size in bp")
        parser.add_argument("--cluster-distance", type=int,
                            default=DEFAULT_CLUSTER_DISTANCE,
                            help="Maximum distance to collapse pairs of break-end calls")
        if len(sys.argv) <= 1:
            parser.parse_args(["--help"])
            sys.exit(0)
        args = parser.parse_args()
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
        meis_file = args.input
        vcf_path = args.output
        sample_name = args.sample_name
        min_score = args.min_alignment_score
        min_reads = args.min_clipped_reads
        min_pct_length = args.min_alignment_percent_length
        min_pct_identity = args.min_alignment_percent_identity
        alu_size = args.alu_size
        sva_size = args.sva_size
        l1_size = args.l1_size
        cluster_distance = args.cluster_distance

    with open(vcf_path, "w") as vcf_out:
        process_meis_to_vcf(
            meis_file, vcf_out, sample_name,
            min_score=min_score, min_reads=min_reads,
            min_pct_length=min_pct_length, min_pct_identity=min_pct_identity,
            alu_size=alu_size, sva_size=sva_size, l1_size=l1_size,
            cluster_distance=cluster_distance,
        )


if __name__ == "__main__":
    main()
