from pathlib import Path

from VCFDataFrame import VCFDataFrame

import pytest

TEST_DATA_PATH = Path("tests/test_data")


def test_load_type():
    with pytest.raises(TypeError) as excinfo:
        VCFDataFrame.read_vcf(18)
    assert "argument must be a string, path to a VCF File" in str(
        excinfo.value
    )


def test_parse_vcf(parsed_vcf):
    df = VCFDataFrame.read_vcf(str(TEST_DATA_PATH / "TEST.vcf"))
    assert df.equals(parsed_vcf)


def test_panel_vcf(vcf_panel):
    df = VCFDataFrame.read_vcf(str(TEST_DATA_PATH / "TEST.vcf"))
    df_panel = df.panel(str(TEST_DATA_PATH / "GeneList.xlsx"))
    assert df_panel.equals(vcf_panel)
