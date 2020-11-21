"""VCFDataFrame main class."""
from VCFDataFrame.io import _load_vcf

import pandas as pd


class VCFDataFrame(pd.DataFrame):
    """Class to extend pandas dataframe using Variant Calling Format."""

    _metadata = ["name"]

    @property
    def _constructor(self):
        """Construct VCFDataFrame class."""
        return VCFDataFrame

    @classmethod
    def read_vcf(cls, vcf):
        """Read VCF file to Pandas DataFrame.

        Parameters
        ----------
        vcf: str
            Path to vcf file

        Returns
        -------
        vcf_df: VCFDataFrame
            VCF file in DataFrame format
        """
        vcf_df, sample_name = _load_vcf(vcf)
        vcf_df = vcf_df.pipe(VCFDataFrame)
        vcf_df.name = sample_name
        return vcf_df

    def priorize_annotations(self):
        """Priorize annotation with a severity basis."""
        pass

    def treat_alt_alleles(self):
        """Treat cases with two or more alternative alleles."""
        pass

    def aminochange(self):
        """Calculate if there is a change in aminoacids, from HGVS.P."""
        pass

    def map_zigosity(self):
        """Map HET and HOM in ZIGOSIS columns."""
        pass

    def process_ESP6500(self):
        """Calculate frequencies on ESP6500 populations and POLYPHEN Score."""
        pass

    def process_CLINVAR(self):
        """Process and translate CLINVAR data."""
        pass
