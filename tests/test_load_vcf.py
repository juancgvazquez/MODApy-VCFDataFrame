from VCFDataFrame import VCFDataFrame

import pytest


def test_transformer_callable():
    with pytest.raises(TypeError) as excinfo:
        VCFDataFrame.read_vcf(18)
    assert "argument must be a string, path to a VCF File" in str(
        excinfo.value
    )
