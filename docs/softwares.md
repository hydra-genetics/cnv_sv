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

