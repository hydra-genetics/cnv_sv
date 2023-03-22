from os.path import abspath

bam_path = abspath(snakemake.input.bam)

with open(snakemake.output.manifest, 'w') as outfile:
    print(bam_path, file=outfile)
