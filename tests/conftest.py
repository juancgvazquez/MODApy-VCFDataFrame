from pathlib import Path

import joblib

import pytest


TEST_DATA_PATH = Path("tests/test_data/")


@pytest.fixture(scope="session")
def parsed_vcf():
    return joblib.load(TEST_DATA_PATH / "parsed_vcf.joblib")


@pytest.fixture(scope="session")
def vcf_panel():
    return joblib.load(TEST_DATA_PATH / "vcf_panel.joblib")
