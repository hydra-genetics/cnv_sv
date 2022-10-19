from pysam import VariantFile
import json

def get_vcf_locus_list(vcf, catalog_list):
    locus_list = []
    str0_loci = []
    for rec in vcf:
        locus_id = rec.info['REPID']
        if rec.alts is not None and rec.alts[0] == '<STR0>':
            str0_loci.append(locus_id.split('_')[0])
            continue
        if rec.filter.get('PASS') is not None:
            if locus_id not in catalog_list:
                locus_id = locus_id.split('_')[0] # extract the locus id from rep id
                if locus_id in catalog_list and locus_id not in locus_list:
                    locus_list.append(locus_id)
            else:
                if locus_id not in locus_list:
                    locus_list.append(locus_id)

    # filter out any loci that have str0 as alt allele, as reviewer throws an error for them
    locus_list = [loc for loc in locus_list if loc not in str0_loci]
    print(locus_list)
    return locus_list

def write_locus_file(output, locus_list):
    locus_str = ','.join(locus_list)
    with open(output, 'w') as outfile:
        print(locus_str, file=outfile)

def get_catalog_loci(catalog):
    with open(catalog, 'r') as cc:
        catalog = json.load(cc)
        locus_list = [j['LocusId'] for j in catalog]
    return locus_list

if __name__ == '__main__':
    catalog_list = get_catalog_loci(snakemake.input.cat)
    vcf = VariantFile(snakemake.input.vcf)
    locus_list = get_vcf_locus_list(vcf, catalog_list)
    write_locus_file(snakemake.output.txt, locus_list)
