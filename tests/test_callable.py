import pytest

from VCFDataFrame import VCFDataFrame


def test_transformer_callable():
    with pytest.raises(TypeError) as excinfo:
        VCFDataFrame.from_vcf(18)
    assert "argument must be a string, path to a VCF File" in str(
        excinfo.value
    )
