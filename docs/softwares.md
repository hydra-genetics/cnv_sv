# Softwares used in the cnv_sv module

## [AutoMap](https://github.com/mquinodo/AutoMap)
Tool to find regions of homozygosity (ROHs) from sequencing data of human samples. Used in analysis of germline samples. Takes a `.vcf` file as input and generates a list of ROHs.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__automap__automap#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__automap__automap#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__automap#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__automap#

---

## [CNVkit batch](https://github.com/etal/cnvkit)
CNVkit calls copy number variation in cancer samples. The program uses a panel of normal to correct for biases and adjust calls based on estimated tumor content (external data). This rule produces copy number segments that are not adjusted for tumor content.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__cnvkit__cnvkit_batch#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__cnvkit__cnvkit_batch#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__cnvkit_batch#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__cnvkit_batch#

---

## [CNVkit call](https://github.com/etal/cnvkit)
CNVkit calls copy number variation in cancer samples. The program uses a panel of normal to correct for biases and adjust calls based on estimated tumor content (external data). This rule produces copy number segments that are adjusted for tumor content.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__cnvkit__cnvkit_call#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__cnvkit__cnvkit_call#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__cnvkit_call#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__cnvkit_call#

---

## [CNVkit diagram](https://github.com/etal/cnvkit)
CNVkit calls copy number variation in cancer samples. The program uses a panel of normal to correct for biases and adjust calls based on estimated tumor content (external data). This rule produces a overview plot of each chromosome with copy number aberrations color coded based on type of variation.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__cnvkit__cnvkit_diagram#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__cnvkit__cnvkit_diagram#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__cnvkit_diagram#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__cnvkit_diagram#

---

## [CNVkit scatter](https://github.com/etal/cnvkit)
CNVkit calls copy number variation in cancer samples. The program uses a panel of normal to correct for biases and adjust calls based on estimated tumor content (external data). This rule produces a copy number plot over the entire genome.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__cnvkit__cnvkit_scatter#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__cnvkit__cnvkit_scatter#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__cnvkit_scatter#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__cnvkit_scatter#

---

## [CNVkit vcf](https://github.com/etal/cnvkit)
CNVkit calls copy number variation in cancer samples. The program uses a panel of normal to correct for biases and adjust calls based on estimated tumor content (external data). This rule exports the called segments into a 'vcf' file used in downstream analysis.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__cnvkit__cnvkit_vcf#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__cnvkit__cnvkit_vcf#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__cnvkit_vcf#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__cnvkit_vcf#

---

## [CNVkit export seg](https://github.com/etal/cnvkit)
CNVkit calls copy number variation in cancer samples. The program uses a panel of normal to correct for biases and adjust calls based on estimated tumor content (external data). This rule exports the called segments into a 'seg' file used in downstream analysis.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__cnvkit__cnvkit_export_seg#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__cnvkit__cnvkit_export_seg#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__cnvkit_export_seg#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__cnvkit_export_seg#

---

## [CNVpytor readdepth](https://github.com/abyzovlab/CNVpytor)
CNVpytor calls copy number variation in WGS germline samples. This rule calculates read depths and creates a binary 'pytor' file containing copy number data.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__cnvpytor__cnvpytor_readdepth#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__cnvpytor__cnvpytor_readdepth#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__cnvpytor_readdepth#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__cnvpytor_readdepth#

---

## [CNVpytor filter](https://github.com/abyzovlab/CNVpytor)
CNVpytor calls copy number variation in WGS germline samples. This rule creates filtered and unfiltered 'vcf' files with called CNVs.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__cnvpytor__cnvpytor_filter#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__cnvpytor__cnvpytor_filter#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__cnvpytor_filter#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__cnvpytor_filter#

---

## [ExomeDepth call](https://github.com/abyzovlab/CNVpytor)
ExomeDepth is a R package designed to detect inherited copy number variants (CNVs) using high throughput DNA sequence data (WES or panles). This rule calls the CNVs.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__exomedepth__exomedepth_call#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__exomedepth__exomedepth_call#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__exomedepth_call#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__exomedepth_call#

---

## [ExomeDepth sex](https://github.com/abyzovlab/CNVpytor)
ExomeDepth is a R package designed to detect inherited copy number variants (CNVs) using high throughput DNA sequence data (WES or panles). This rule only copies the peddy sex file to a copy used by ExomeDepth.

### :snake: Rule

SNAKEMAKE_RULE_SOURCE__exomedepth__exomedepth_sex#

#### :left_right_arrow: input / output files

SNAKEMAKE_RULE_TABLE__exomedepth__exomedepth_sex#

### :wrench: Configuration

#### Software settings (`config.yaml`)

CONFIGSCHEMA__exomedepth_sex#

#### Resources settings (`resources.yaml`)

RESOURCESSCHEMA__exomedepth_sex#

---

## [ExpansionHunter](https://github.com/Illumina/ExpansionHunter)
Expansion Hunter aims to estimate sizes of selected repeats by performing a targeted search through a BAM/CRAM file for reads that span, flank, and are fully contained in each repeat.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__expansionhunter__expansionhunter#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__expansionhunter__expansionhunter#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__expansionhunter#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__expansionhunter#

---

## [gatk_collect_read_counts](https://gatk.broadinstitute.org/hc/en-us/articles/360037592671-CollectReadCounts)
A collection of rules for GATK CNV calling in somatic samples using a panel of normal following GATK [guide](https://gatk.broadinstitute.org/hc/en-us/articles/360035531092--How-to-part-I-Sensitively-detect-copy-ratio-alterations-and-allelic-segments). This rule goes through the `bam` file of the sample and collects coverage statistics.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__gatk__gatk_collect_read_counts#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__gatk__gatk_collect_read_counts#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__gatk_collect_read_counts#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__gatk_collect_read_counts#

---

## [gatk_collect_allelic_counts](https://gatk.broadinstitute.org/hc/en-us/articles/360037594071-CollectAllelicCounts)
A collection of rules for GATK CNV calling in somatic samples using a panel of normal following GATK [guide](https://gatk.broadinstitute.org/hc/en-us/articles/360035531092--How-to-part-I-Sensitively-detect-copy-ratio-alterations-and-allelic-segments). This rule collects allele counts for germline SNPs specified in input.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__gatk__gatk_collect_allelic_counts#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__gatk__gatk_collect_allelic_counts#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__gatk_collect_allelic_counts#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__gatk_collect_allelic_counts#

---

## [gatk_denoise_read_counts](https://gatk.broadinstitute.org/hc/en-us/articles/360040508731-DenoiseReadCounts)
A collection of rules for GATK CNV calling in somatic samples using a panel of normal following GATK [guide](https://gatk.broadinstitute.org/hc/en-us/articles/360035531092--How-to-part-I-Sensitively-detect-copy-ratio-alterations-and-allelic-segments). This rule uses a panel of normal to corrects the coverage data collected by the rule gatk_collect_read_counts.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__gatk__gatk_denoise_read_counts#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__gatk__gatk_denoise_read_counts#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__gatk_denoise_read_counts#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__gatk_denoise_read_counts#

---

## [gatk_model_segments](https://gatk.broadinstitute.org/hc/en-us/articles/360037225892-ModelSegments)
A collection of rules for GATK CNV calling in somatic samples using a panel of normal following GATK [guide](https://gatk.broadinstitute.org/hc/en-us/articles/360035531092--How-to-part-I-Sensitively-detect-copy-ratio-alterations-and-allelic-segments). This rule makes the actual copy number segmentation.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__gatk__gatk_model_segments#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__gatk__gatk_model_segments#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__gatk_model_segments#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__gatk_model_segments#

---

## [gatk_call_copy_ratio_segments](https://gatk.broadinstitute.org/hc/en-us/articles/360037593771-CallCopyRatioSegments)
A collection of rules for GATK CNV calling in somatic samples using a panel of normal following GATK [guide](https://gatk.broadinstitute.org/hc/en-us/articles/360035531092--How-to-part-I-Sensitively-detect-copy-ratio-alterations-and-allelic-segments). Calls amplifications and deletion based on a statistical test.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__gatk__gatk_call_copy_ratio_segments#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__gatk__gatk_call_copy_ratio_segments#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__gatk_call_copy_ratio_segments#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__gatk_call_copy_ratio_segments#

---

## [gatk_to_vcf](https://github.com/hydra-genetics/cnv_sv/blob/develop/workflow/scripts/gatk_to_vcf.py)
An in-house script that converts the GATK segmentation into a file compatible with the `vcf` file format.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__gatk__gatk_to_vcf#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__gatk__gatk_to_vcf#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__gatk_to_vcf#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__gatk_to_vcf#

---

## [manta config tn](https://github.com/Illumina/manta)
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. This rule configures a python script that runs Manta to call somatic variants in tumor normal mode. 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__manta__manta_config_tn#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__manta__manta_config_tn#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__manta_config_tn#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__manta_config_tn#

---

## [manta run workflow tn](https://github.com/Illumina/manta)
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. This rule runs Manta to call somatic variants in tumor normal mode. 


### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__manta__manta_run_workflow_tn#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__manta__manta_run_workflow_tn#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__manta_run_workflow_tn#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__manta_run_workflow_tn#

---

## [manta config t](https://github.com/Illumina/manta)
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. This rule configures a python script that runs Manta to call somatic variants in tumor only mode. 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__manta__manta_config_t#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__manta__manta_config_t#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__manta_config_t#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__manta_config_t#

---

## [manta run workflow t](https://github.com/Illumina/manta)
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. This rule runs Manta to call somatic variants in tumor only mode. 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__manta__manta_run_workflow_t#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__manta__manta_run_workflow_t#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__manta_run_workflow_t#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__manta_run_workflow_t#

---

## [manta config n](https://github.com/Illumina/manta)
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. This rule configures a python script that runs Manta to call germline variants in normal samples. 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__manta__manta_config_n#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__manta__manta_config_n#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__manta_config_n#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__manta_config_n#

---

## [manta run workflow n](https://github.com/Illumina/manta)
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. This rule runs Manta to call germline variants in normal samples. 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__manta__manta_run_workflow_n#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__manta__manta_run_workflow_n#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__manta_run_workflow_n#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__manta_run_workflow_n#

---