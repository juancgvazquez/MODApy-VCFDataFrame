"""Supports every IO operation."""
import logging
import os
from collections import OrderedDict

import cyvcf2

import pandas as pd


def _load_vcf(vcf):
    """VCF Parser to a pd.DataFrame.

    Parameters
    ----------
    vcf
        Path to the vcf to parse.

    Returns
    -------
    vcf_df
        DataFrame created from vcf
    name
        Sample Name
    """
    if not isinstance(vcf, str):
        logging.error(f"Received argument was {vcf}.")
        raise TypeError("argument must be a string, path to a VCF File")
    if not vcf.lower().endswith(".vcf"):
        logging.error(f"Received argument was {vcf}.")
        raise TypeError("filepath must end with .vcf")
    if not os.path.exists(vcf):
        logging.error(f"Received argument was {vcf}.")
        raise FileNotFoundError("File not found in vcf path")
    logging.info(f"Parsing VCF File: {vcf}")
    vcf_reads = cyvcf2.Reader(vcf)
    if len(vcf_reads.samples) > 0:
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
    return vcf_df, name, vcf_reads


def _load_panel(panel):
    try:
        paneldf = pd.ExcelFile(panel).parse("GeneList")
        gene_symbol_list = list(paneldf.GeneSymbol.unique())
    except Exception as err:
        raise Exception(
            f"There was an error parsing GeneList.\
                        Error was: {err}"
        )
    return gene_symbol_list
