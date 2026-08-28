import pytest
from types import SimpleNamespace
from hydra_genetics.utils.misc import get_input_haplotagged_bam


# Replicated here since common.smk is not importable as a Python module.
# Must be kept in sync with workflow/rules/common.smk: get_severus_tn_input
def get_severus_tn_input(wildcards, config):
    bam_t, bai_t = get_input_haplotagged_bam(wildcards, config, set_type="T")
    bam_n, bai_n = get_input_haplotagged_bam(wildcards, config, set_type="N")
    return {
        "bam_t": bam_t,
        "bai_t": bai_t,
        "bam_n": bam_n,
        "bai_n": bai_n,
    }


@pytest.fixture
def wildcards():
    return SimpleNamespace(sample="NA12878", type="T")


def test_default_config_returns_whatshap_haplotag_paths(wildcards):
    """With no 'phaser' set, the upstream default path is used."""
    config = {}
    result = get_severus_tn_input(wildcards, config)
    assert result["bam_t"] == "snv_indels/whatshap_haplotag/NA12878_T.haplotagged.bam"
    assert result["bai_t"] == "snv_indels/whatshap_haplotag/NA12878_T.haplotagged.bam.bai"
    assert result["bam_n"] == "snv_indels/whatshap_haplotag/NA12878_N.haplotagged.bam"
    assert result["bai_n"] == "snv_indels/whatshap_haplotag/NA12878_N.haplotagged.bam.bai"


def test_wildcard_type_is_ignored(wildcards):
    """bam_t must always be T and bam_n must always be N, regardless of wildcards.type."""
    wildcards.type = "N"
    config = {}
    result = get_severus_tn_input(wildcards, config)
    assert result["bam_t"] == "snv_indels/whatshap_haplotag/NA12878_T.haplotagged.bam"
    assert result["bam_n"] == "snv_indels/whatshap_haplotag/NA12878_N.haplotagged.bam"


def test_phaser_whatshap_is_respected(wildcards):
    config = {"phaser": "whatshap"}
    result = get_severus_tn_input(wildcards, config)
    assert result["bam_t"] == "snv_indels/whatshap_haplotag/NA12878_T.haplotagged.bam"
    assert result["bam_n"] == "snv_indels/whatshap_haplotag/NA12878_N.haplotagged.bam"


def test_phaser_hiphase_is_respected(wildcards):
    """This is the phaser used by .tests/integration/config_pacbio.yaml."""
    config = {"phaser": "hiphase"}
    result = get_severus_tn_input(wildcards, config)
    assert result["bam_t"] == "snv_indels/hiphase/NA12878_T.haplotagged.bam"
    assert result["bai_t"] == "snv_indels/hiphase/NA12878_T.haplotagged.bam.bai"
    assert result["bam_n"] == "snv_indels/hiphase/NA12878_N.haplotagged.bam"
    assert result["bai_n"] == "snv_indels/hiphase/NA12878_N.haplotagged.bam.bai"


def test_unknown_phaser_falls_back_to_snv_indels_prefix(wildcards):
    """An unmapped phaser name is used verbatim under snv_indels/."""
    config = {"phaser": "longphase"}
    result = get_severus_tn_input(wildcards, config)
    assert result["bam_t"] == "snv_indels/longphase/NA12878_T.haplotagged.bam"
    assert result["bam_n"] == "snv_indels/longphase/NA12878_N.haplotagged.bam"


def test_returns_dict_with_expected_keys(wildcards):
    result = get_severus_tn_input(wildcards, {})
    assert set(result.keys()) == {"bam_t", "bai_t", "bam_n", "bai_n"}
