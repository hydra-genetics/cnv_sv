import pytest
from types import SimpleNamespace
from hydra_genetics.utils.misc import get_input_haplotagged_bam


# Replicated here since common.smk is not importable as a Python module.
# Must be kept in sync with workflow/rules/common.smk: get_severus_tn_input
def get_severus_tn_input(wildcards, config):
    wc_t = SimpleNamespace(sample=wildcards.sample, type="T")
    wc_n = SimpleNamespace(sample=wildcards.sample, type="N")
    bam_t, bai_t = get_input_haplotagged_bam(wc_t, config)
    bam_n, bai_n = get_input_haplotagged_bam(wc_n, config)
    return {
        "bam_t": bam_t,
        "bai_t": bai_t,
        "bam_n": bam_n,
        "bai_n": bai_n,
    }


@pytest.fixture
def wildcards():
    return SimpleNamespace(sample="NA12878", type="T")


def test_default_config_returns_merge_bam_paths(wildcards):
    config = {}
    result = get_severus_tn_input(wildcards, config)
    assert result["bam_t"] == "alignment/samtools_merge_bam/NA12878_T.bam"
    assert result["bai_t"] == "alignment/samtools_merge_bam/NA12878_T.bam.bai"
    assert result["bam_n"] == "alignment/samtools_merge_bam/NA12878_N.bam"
    assert result["bai_n"] == "alignment/samtools_merge_bam/NA12878_N.bam.bai"


def test_wildcard_type_is_ignored(wildcards):
    """bam_t must always be T and bam_n must always be N, regardless of wildcards.type."""
    wildcards.type = "N"
    config = {}
    result = get_severus_tn_input(wildcards, config)
    assert result["bam_t"] == "alignment/samtools_merge_bam/NA12878_T.bam"
    assert result["bam_n"] == "alignment/samtools_merge_bam/NA12878_N.bam"


def test_haplotag_path_is_respected(wildcards):
    config = {"haplotag_path": "alignment/longphase_haplotag"}
    result = get_severus_tn_input(wildcards, config)
    assert result["bam_t"] == "alignment/longphase_haplotag/NA12878_T.bam"
    assert result["bai_t"] == "alignment/longphase_haplotag/NA12878_T.bam.bai"
    assert result["bam_n"] == "alignment/longphase_haplotag/NA12878_N.bam"
    assert result["bai_n"] == "alignment/longphase_haplotag/NA12878_N.bam.bai"


def test_haplotag_suffix_is_respected(wildcards):
    config = {"haplotag_suffix": "haplotagged"}
    result = get_severus_tn_input(wildcards, config)
    assert result["bam_t"] == "alignment/samtools_merge_bam/NA12878_T.haplotagged.bam"
    assert result["bai_t"] == "alignment/samtools_merge_bam/NA12878_T.haplotagged.bam.bai"
    assert result["bam_n"] == "alignment/samtools_merge_bam/NA12878_N.haplotagged.bam"
    assert result["bai_n"] == "alignment/samtools_merge_bam/NA12878_N.haplotagged.bam.bai"


def test_haplotag_path_and_suffix_combined(wildcards):
    config = {"haplotag_path": "alignment/longphase_haplotag", "haplotag_suffix": "haplotagged"}
    result = get_severus_tn_input(wildcards, config)
    assert result["bam_t"] == "alignment/longphase_haplotag/NA12878_T.haplotagged.bam"
    assert result["bai_t"] == "alignment/longphase_haplotag/NA12878_T.haplotagged.bam.bai"
    assert result["bam_n"] == "alignment/longphase_haplotag/NA12878_N.haplotagged.bam"
    assert result["bai_n"] == "alignment/longphase_haplotag/NA12878_N.haplotagged.bam.bai"


def test_returns_dict_with_expected_keys(wildcards):
    result = get_severus_tn_input(wildcards, {})
    assert set(result.keys()) == {"bam_t", "bai_t", "bam_n", "bai_n"}
