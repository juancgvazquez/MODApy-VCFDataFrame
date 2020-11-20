import logging
import os
from collections import OrderedDict

import cyvcf2
import pandas as pd


class VCFDataFrame(pd.DataFrame):
    """Class to extend pandas dataframe using Variant Calling Format"""

    _metadata = ["name"]

    @property
    def _constructor(self):
        """Method to mantain class type"""
        return VCFDataFrame

    @classmethod
    def from_vcf(cls, vcf):
        """
        Method that creates a ParsedVCF1 (a DataFrame) from a vcf file
        Parameters
        ----------
        vcf
            Path to the vcf to parse.
        """
        if not isinstance(vcf, str):
            raise TypeError("argument must be a string, path to a VCF File")
        if not vcf.lower().endswith(".vcf"):
            raise TypeError("filepath must end with .vcf")
        if not os.exists(vcf):
            raise FileNotFoundError("File not found in vcf path")
        logging.info(f"Parsing VCF File: {vcf}")
        vcf_reads = cyvcf2.Reader(vcf)
        if vcf_reads.samples > 0:
            name = vcf_reads.samples[0]
        else:
            name = vcf.split("/")[-1]
        variantsDict = OrderedDict()
        for variant in vcf_reads:
            variantsDict[
                variant.CHROM
                + "+"
                + str(variant.POS)
                + "+"
                + variant.REF
                + "+"
                + ",".join(variant.ALT)
            ] = {
                "ID": variant.ID,
                "QUAL": variant.QUAL,
                "FILTER": variant.FILTER,
            }
            variantsDict[
                variant.CHROM
                + "+"
                + str(variant.POS)
                + "+"
                + variant.REF
                + "+"
                + ",".join(variant.ALT)
            ].update({k: v for (k, v) in variant.INFO})

        vcf_df = pd.DataFrame.from_dict(variantsDict, orient="index")
        del variantsDict
        vcf_df.index = vcf_df.index.str.split("+", expand=True)
        vcf_df.index.names = ["CHROM", "POS", "REF", "ALT"]
        vcf_df.reset_index(inplace=True)
        vcf_df = vcf_df.pipe(VCFDataFrame)
        vcf_df.name = name
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
