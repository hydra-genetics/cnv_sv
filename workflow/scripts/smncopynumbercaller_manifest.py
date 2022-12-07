from os.path import abspath

bam_path = abspath(snakemake.input.bam)

print(bam_path)
with open(snakemake.output.manifest, 'w') as outfile:
    print(bam_path, file=outfile)
