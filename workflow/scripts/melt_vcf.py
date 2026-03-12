#!/usr/bin/env python3
"""Fix MELT VCF:
1. Set GL Number=G in FORMAT header.
2. Replace ME-specific ALT headers (INS:ME:*) with a single ##ALT=<ID=INS,...> header.
3. Normalize record ALT to <INS>.
4. Normalize record SVTYPE to INS.
"""

__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2026, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"

import gzip
import re


def fix_melt_vcf(input_vcf, output_vcf):
    alt_header_written = False

    opener = gzip.open if input_vcf.endswith(".gz") else open

    with opener(input_vcf, "rt") as fin, open(output_vcf, "w") as fout:
        for line in fin:
            if line.startswith("##FORMAT=<ID=GL,"):
                fout.write(re.sub(r"Number=\d+", "Number=G", line))
            elif line.startswith("##ALT=<ID=INS:ME:"):
                if not alt_header_written:
                    fout.write('##ALT=<ID=INS,Description="Insertion">\n')
                    alt_header_written = True
            elif line.startswith("#"):
                fout.write(line)
            else:
                fields = line.split("\t")
                fields[4] = "<INS>"
                fields[7] = re.sub(r"SVTYPE=[^;]+", "SVTYPE=INS", fields[7])
                fout.write("\t".join(fields))


if __name__ == "__main__":
    fix_melt_vcf(snakemake.input.vcf, snakemake.output.vcf)
