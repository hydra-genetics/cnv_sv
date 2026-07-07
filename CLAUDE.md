# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the **hydra-genetics/cnv_sv** Snakemake workflow module for calling Copy Number Variants (CNVs) and Structural Variants (SVs). It supports both short-read (Illumina) and long-read (PacBio/Oxford Nanopore) sequencing platforms, and is designed to be used as a module within the broader hydra-genetics ecosystem.

## Commands

### Workflow Validation (from `.tests/integration/`)
```bash
# Dry-run for short-read config
snakemake -n -s ../../workflow/Snakefile --configfile config.yaml

# Dry-run for long-read/PacBio config
snakemake -n -s ../../workflow/Snakefile --configfile config_pacbio.yaml

# Lint the Snakefile
snakemake --lint -n -s ../../workflow/Snakefile --configfile config.yaml
```

### Python Tests
```bash
# Run all script unit tests
pytest workflow/scripts/

# Run a single test file
pytest workflow/scripts/scramble_vcf_test.py
```

### Code Style
```bash
# Check Snakemake formatting (line length 130, skip string normalization)
snakefmt workflow/rules/

# Check Python style (max line length 130)
pycodestyle --max-line-length=130 workflow/scripts/
```

## Architecture

### Entry Point
`workflow/Snakefile` is the main entry point. It includes all rule files and defines:
- Wildcard constraints for sample/type/flowcell/lane/barcode patterns
- Rule ordering (e.g., `ruleorder: trgt_genotype > bgzip`)
- The `rule all` that compiles expected outputs via `compile_output_list()`

### Key Rule Files
- `workflow/rules/common.smk` — Shared helper functions, wildcard constraints, and the `compile_output_list()` function that drives all output generation. Read this first when debugging missing outputs.
- Tool-specific rules are in `workflow/rules/*.smk` (one file per tool: `cnvkit.smk`, `manta.smk`, `severus.smk`, etc.)

### Configuration
- `workflow/schemas/config.schema.yaml` — Defines all valid config keys and their types. Check here before adding new parameters.
- `.tests/integration/config.yaml` — Short-read integration test config (reference for how real configs look)
- `.tests/integration/config_pacbio.yaml` — Long-read integration test config

### Data Flow
1. **Inputs**: BAM files from an upstream `alignment` module (hydra-genetics convention)
2. **CNV/SV calling**: Multiple independent callers run in parallel
3. **Merging**: Results merged via `svdb_merge` (see `workflow/rules/svdb.smk`)
4. **Tumor analysis**: Optional tumor purity estimation via PureCN (`workflow/rules/purecn.smk`)

### Platform Awareness
Rules check `units.platform` to distinguish short-read vs. long-read:
- Short-read: platform NOT in `["ONT", "PACBIO"]`
- Long-read: platform in `["ONT", "PACBIO"]`
- Short-read sample indexing: `[sample, type, flowcell, lane, barcode]`
- Long-read sample indexing: `[sample, type, processing_unit, barcode]`

### Rule Structure Convention
Every rule follows this pattern:
```python
rule rule_name:
    input:    # BAM files or upstream rule outputs
    output:   # Usually temp() for intermediates
    params:   # From config.yaml via config.get()
    log:
    benchmark:
    threads:
    resources:
    container: # Docker/Singularity image
    message:
    shell:
```

### Python Scripts
`workflow/scripts/` contains data transformation scripts. Unit tests live alongside scripts with a `_test.py` suffix (e.g., `scramble_vcf.py` + `scramble_vcf_test.py`).

## Style Guidelines
- Max line length: **130 characters** (both Python and Snakemake)
- Snakemake files formatted with `snakefmt` (skip string normalization)
- Python style enforced with `pycodestyle`

## Key Helper Functions (common.smk)
- `compile_output_list()` — Generates all expected output files based on samples/units data
- `get_tc()` — Gets tumor content from pathology data or computed sources
- `get_karyotype()` — Maps sex to karyotype string
- `get_vcfs_for_svdb_merge()` — Conditionally aggregates VCFs for merging
- `get_purecn_inputs()` — Resolves dynamic inputs based on segmentation method

Other helper functions (like, `get_input_aligned_bam()`, `get_input_haplotagged_bam()`, etc.) are imported from the main hydra-genetics module (https://github.com/hydra-genetics/hydra-genetics/blob/develop/hydra_genetics/utils/misc.py) at runtime.
