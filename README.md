<p align="center">
<h1><a href="https://hydra-genetics-cnv-sv.readthedocs.io">https://hydra-genetics-cnv-sv.readthedocs.io</a><h1>
</p>

# <img src="images/hydragenetics.png" width=40 /> hydra-genetics/cnv_sv

Snakemake module containing steps to call copy number variants and structural variants

![Lint](https://github.com/hydra-genetics/cnv_sv/actions/workflows/lint.yaml/badge.svg?branch=master)
![Snakefmt](https://github.com/hydra-genetics/cnv_sv/actions/workflows/snakefmt.yaml/badge.svg?branch=master)
![snakemake dry run](https://github.com/hydra-genetics/cnv_sv/actions/workflows/snakemake-dry-run.yaml/badge.svg?branch=master)
![integration test](https://github.com/hydra-genetics/cnv_sv/actions/workflows/integration.yaml/badge.svg?branch=master)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/gpl-3.0.html)

## :speech_balloon: Introduction

The module contain rules used to call CNVs (copy number variants) and SV (structural variants), also a rule used to
merge CNV/SV

## :heavy_exclamation_mark: Dependencies

In order to use this module, the following dependencies are required:

[![hydra-genetics](https://img.shields.io/badge/hydragenetics-v0.15.0-blue)](https://github.com/hydra-genetics/)
[![pandas](https://img.shields.io/badge/pandas-1.3.1-blue)](https://pandas.pydata.org/)
[![python](https://img.shields.io/badge/python-3.8-blue)](https://www.python.org/)
[![snakemake](https://img.shields.io/badge/snakemake-7.8.0-blue)](https://snakemake.readthedocs.io/en/stable/)
[![singularity](https://img.shields.io/badge/singularity-3.0.0-blue)](https://sylabs.io/docs/)

## :school_satchel: Preparations

### Sample and unit data short read

Input data should be added to [`samples.tsv`](https://github.com/hydra-genetics/prealignment/blob/develop/config/samples.tsv)
and [`units.tsv`](https://github.com/hydra-genetics/prealignment/blob/develop/config/units.tsv).
The following information need to be added to these files:

| Column Id | Description |
| --- | --- |
| **`samples.tsv`** |
| sample | unique sample/patient id, one per row |
| tumor_content | ratio of tumor cells to total cells |
| **`units.tsv`** |
| sample | same sample/patient id as in `samples.tsv` |
| type | data type identifier (one letter), can be one of **T**umor, **N**ormal, **R**NA |
| platform | type of sequencing platform, e.g. `NovaSeq` |
| machine | specific machine id, e.g. NovaSeq instruments have `@Axxxxx` |
| flowcell | identifer of flowcell used |
| lane | flowcell lane number |
| barcode | sequence library barcode/index, connect forward and reverse indices by `+`, e.g. `ATGC+ATGC` |
| fastq1/2 | absolute path to forward and reverse reads |
| adapter | adapter sequences to be trimmed, separated by comma |


### Sample and unit data long read

Input data should be added to [`samples.tsv`](https://github.com/hydra-genetics/cnv_sv/blob/develop/.tests/integration/samples_pacbio.tsv)
and [`units.tsv`](https://github.com/hydra-genetics/cnv_sv/blob/develop/.tests/integration/units_pacbio.tsv).
The following information need to be added to these files:

| Column Id | Description |
| --- | --- |
| **`samples.tsv`** |
| sample | unique sample/patient id, one per row |
| tumor_content | ratio of tumor cells to total cells |
| **`units.tsv`** |
| sample | same sample/patient id as in `samples.tsv` |
| type | data type identifier (one letter), can be one of **T**umor, **N**ormal, **R**NA |
| platform | type of sequencing platform, e.g. `Pacbio` |
| machine | specific machine id, e.g. Pacbio instrument e.g. `REVIO` |
| processing_unit | contains the PacBio movie name or flowcell id for ONT |
| barcode | sequence library barcode/index, connect forward and reverse indices by `-`, e.g. `ATGC-ATGC` |
| methylation | Methylation read tags are found in the bam file (Yes or No) |
| bam | absolute path to the unmapped BAM file |


### Reference data

A reference .fasta-file should be specified in config.yaml in the section reference and fasta. In addition,
the file should be indexed using samtools faidx and the path of the resulting file added to the stanza fai.
A bed file containing the covered regions shall be added to design_bed.

### Panel of normals (PoN) and read count

PoN must be configured to be able to run CNVkit and gatk CNV, workflows for this can be found at [hydra-genetics/references](
https://github.com/hydra-genetics/references). Instructions for the tools can be found at.
 * [CNVkit](https://cnvkit.readthedocs.io/en/stable/quickstart.html?highlight=panel%20of%20normal#build-a-reference-from-normal-samples-and-infer-tumor-copy-ratios)
 * [GATK CNV](https://gatk.broadinstitute.org/hc/en-us/articles/360035531092--How-to-part-I-Sensitively-detect-copy-ratio-alterations-and-allelic-segments#2)

For exomdepth read count data must be generated. more information can be found at [ExomeDepth repo](https://github.com/vplagnol/ExomeDepth/blob/80da0cb76d6a9a0ad4c422ea5a9ff3b82f9f6279/vignette/vignette.Rnw#L114)

### VCF used to calculate variant allele frequencies
Could be a gnomeAD vcf file filtered on population allele frequences above 0.001

### Purecn
Purecn is used to estimate tumor purity from the data. Purecn can be run with different kinds of segmentation methods in conjunction with different variant files as input, see further purecn in the [schemas](https://github.com/hydra-genetics/cnv_sv/blob/develop/workflow/schemas/config.schema.yaml). 

## :white_check_mark: Testing

The workflow repository contains a small test dataset `.tests/integration` which can be run like so:

```bash
$ cd .tests/integration
$ snakemake -s ../../Snakefile -j1 --configfile config.yaml --use-singularity
```

## :rocket: Usage

To use this module in your workflow, follow the description in the
[snakemake docs](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).
Add the module to your `Snakefile` like so:

```bash
module cns_sv:
    snakefile:
        github(
            "hydra-genetics/cnv_sv",
            path="workflow/Snakefile",
            tag="v0.1.0",
        )
    config:
        config


use rule * from cnv_sv as cnv_sv_*
```

### Compatibility

Latest:
 - alignment:v0.4.0

 See [COMPATIBLITY.md](../master/COMPATIBLITY.md) file for a complete list of module compatibility.


### Disclaimer

Since it is not possible to create integration test without a full dataset purecn will not be subjected to integration testing 
and we can not guarantee that it will work.

Since it is not possible to create integration test without a large dataset SMNCopyNumberCaller will not be subjected to integration
testing and we can not guarantee that it will work
