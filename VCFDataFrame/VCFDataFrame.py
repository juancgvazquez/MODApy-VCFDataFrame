"""VCFDataFrame main class."""
from VCFDataFrame.io import _load_panel, _load_vcf

import numpy as np

import pandas as pd


class VCFDataFrame(pd.DataFrame):
    """Class to extend pandas dataframe using Variant Calling Format."""

    _metadata = ["name"]

    @property
    def _constructor(self):
        """Construct VCFDataFrame class."""
        return VCFDataFrame

    @classmethod
    def read_vcf(
        cls,
        vcf,
        priorize_ann=True,
        aminochange=True,
        zigosity=True,
        parse_ESP=True,
        parse_CLINVAR=True,
        round_nums=6,
    ):
        """Read VCF file to Pandas DataFrame.

        Parameters
        ----------
        vcf: str
            Path to vcf file
        priorize_ann: bool, dict
            Priorize annotations, receives bool or priorization dict
        aminochange: bool
            Calculate AminoChange column
        zigosity: bool
            Parse zigosity and map to HOM or HET
        parse_ESP: bool
            Parse frequencies in ESP6500 annotations
        parse_CLINVAR: bool
            Parse and translate CLINVAR annotation codes
        round_nums: int or None
            Round nums to N decimals. None value will stop round method

        Returns
        -------
        vcf_df: VCFDataFrame
            VCF file in DataFrame format
        """
        vcf_df, sample_name, pVCF = _load_vcf(vcf)
        vcf_df = vcf_df.pipe(VCFDataFrame)
        vcf_df = vcf_df._treat_alt_alleles()
        if priorize_ann is True:
            vcf_df = vcf_df._priorize_annotations(pVCF)
        elif isinstance(priorize_ann, dict):
            vcf_df = vcf_df._priorize_annotations(pVCF, priorize_ann)
        vcf_df.columns = vcf_df.columns.str.upper()
        if aminochange is True:
            vcf_df._calc_aminochange()
        if zigosity is True:
            vcf_df = vcf_df._map_zigosity()
        if parse_ESP is True:
            vcf_df = vcf_df._process_ESP6500()
        if parse_CLINVAR is True:
            vcf_df = vcf_df._process_CLINVAR()
        if round_nums is not None:
            vcf_df = vcf_df._round_num_cols(pVCF, round_nums)
        vcf_df.replace(["nan", "", np.nan], ".", inplace=True)
        vcf_df.replace([[None], "."], inplace=True, regex=True)
        vcf_df = vcf_df.astype("str")
        vcf_df["POS"] = vcf_df["POS"].astype(int)
        vcf_df = vcf_df.pipe(VCFDataFrame)
        vcf_df.name = sample_name
        return vcf_df

    def _round_num_cols(self, pVCF, round_nums):
        numcols = list()
        for x in pVCF.header_iter():
            if x.type == "INFO":
                if x["Type"] in ["Float", "Integer"]:
                    numcols.append(x["ID"])
        additional_cols = [
            i
            for i in self.columns
            if i in ["ESP6500_MAF_EA", "ESP6500_MAF_AA", "ESP6500_MAF_ALL"]
        ]
        numcols += additional_cols
        numcols = list(
            set(
                [
                    x.upper()
                    for x in numcols
                    for y in self.columns
                    if x.upper() == y
                ]
            )
        )
        self[numcols] = self[numcols].apply(
            pd.to_numeric, errors="coerce", axis=1
        )
        self = self.round(round_nums)
        return self

    def _priorize_annotations(self, pVCF, priority_dict=None):
        """Priorize annotation with a severity basis."""
        if "ANN" in self.columns:
            anndf = self["ANN"]
            annhead = pVCF.get_header_type("ANN")["Description"].strip(
                '"Functional annotations: \'"'
            )
            annheaderlist = [x.strip() for x in annhead.split("|")]
            anndf = anndf.str.split(",", expand=True).stack()
            anndf = anndf.str.split("|", expand=True)
            anndf.columns = annheaderlist
            self.drop(columns="ANN", inplace=True)
            anndf.index = anndf.index.droplevel(1)
            self = self.join(anndf, how="inner")
            del anndf
            del annhead
            del annheaderlist
            if priority_dict is None:
                IMPACT_SEVERITY = {
                    "exon_loss_variant": 1,
                    "frameshift_variant": 2,
                    "stop_gained": 3,
                    "stop_lost": 4,
                    "start_lost": 5,
                    "splice_acceptor_variant": 6,
                    "splice_donor_variant": 7,
                    "disruptive_inframe_deletion": 8,
                    "inframe_insertion": 9,
                    "disruptive_inframe_insertion": 10,
                    "inframe_deletion": 11,
                    "missense_variant": 12,
                    "splice_region_variant": 13,
                    "stop_retained_variant": 14,
                    "initiator_codon_variant": 15,
                    "synonymous_variant": 16,
                    "start_retained": 17,
                    "coding_sequence_variant": 18,
                    "5_prime_UTR_variant": 19,
                    "3_prime_UTR_variant": 20,
                    "5_prime_UTR_premature_start_codon_gain_variant": 21,
                    "intron_variant": 22,
                    "non_coding_exon_variant": 23,
                    "upstream_gene_variant": 24,
                    "downstream_gene_variant": 25,
                    "TF_binding_site_variant": 26,
                    "regulatory_region_variant": 27,
                    "intergenic_region": 28,
                    "transcript": 29,
                }
            else:
                IMPACT_SEVERITY = priority_dict
            self["sorter"] = (
                self["Annotation"]
                .str.split("&")
                .str[0]
                .replace(IMPACT_SEVERITY)
            )
            self.loc[self["HGVS.c"].str.contains("null"), "HGVS.c"] = None
            self["sorter2"] = [
                x[0] == x[1] for x in zip(self["ALT"], self["Allele"])
            ]
            self = self.sort_values(
                by=["CHROM", "POS", "sorter2", "sorter"],
                ascending=[True, True, False, True],
            ).drop_duplicates(["CHROM", "POS", "REF", "ALT"])
            self.drop(columns=["sorter", "sorter2"], inplace=True)
            del IMPACT_SEVERITY
            return self

    def _treat_alt_alleles(self):
        """Treat cases with two or more alternative alleles."""
        splitdf = self.loc[self["ALT"].str.contains(",")].copy()
        if len(splitdf) > 0:
            ALT = (
                splitdf["ALT"]
                .astype(str)
                .str.split(",", n=1, expand=True)
                .stack()
                .rename("ALT")
            )
            ALT.index = ALT.index.droplevel(-1)
            ALT = ALT.to_frame()
            splitdf = splitdf.join(ALT, lsuffix="_x", rsuffix="_y")
            del ALT
            splitdf["ALT"] = splitdf["ALT_y"].combine_first(splitdf["ALT_x"])
            splitdf.drop(columns=["ALT_y", "ALT_x"], inplace=True)
            splitdf.reset_index(inplace=True)
            splitdf.drop(columns="index", inplace=True)
        odd = splitdf.iloc[::2].copy()
        even = splitdf.iloc[1::2].copy()
        splitlist = [
            "ID",
            "AC",
            "AF",
            "SAMPLES_AF",
            "MLEAC",
            "MLEAF",
            "VARTYPE",
            "dbSNPBuildID",
        ]
        splitlist = [x for x in splitlist if x in self.columns]
        splitlist += [
            x for x in self.columns if x.startswith(("1000", "CLINVAR"))
        ]
        for col in splitlist:
            odd[col] = odd[col].astype(str).str.split(",", n=1).str[0]
            even[col] = even[col].apply(
                lambda x: x
                if len(str(x).split(",")) <= 1
                else str(x).split(",", maxsplit=1)[1]
            )
        splitdf = (
            pd.concat([odd, even])
            .sort_index()
            .replace(to_replace=[r"\(", r"\)"], value="", regex=True)
        )
        del odd, even
        splitdf = splitdf[["CHROM", "POS", "REF", "ALT"] + splitlist]
        self = self.merge(splitdf, on=["CHROM", "POS", "REF"], how="left")
        splitlist.append("ALT")
        xlist = [x + "_x" for x in splitlist]
        ylist = [y + "_y" for y in splitlist]
        del splitdf  # ya no uso más splitdf así que lo borro
        for col in splitlist:
            self[col] = self[col + "_y"].combine_first(self[col + "_x"])
        del splitlist  # ya no uso más splitlist
        self.drop(columns=xlist + ylist, inplace=True)
        del xlist, ylist  # ya no uso más esto.
        self["POS"] = self["POS"].astype(int)
        return self

    def _calc_aminochange(self):
        if "HGVS.P" in self.columns:
            self["AMINOCHANGE"] = self["HGVS.P"].apply(self._aminochange)

    @classmethod
    def _aminochange(cls, value):
        """Calculate if there is a change in aminoacids, from HGVS.P."""
        try:
            value = value.replace("p.", "")
            if value[:3] != value[-3:]:
                return "CHANGE"
            else:
                return "."
        except TypeError:
            return "."

    @classmethod
    def _divide(cls, x, y):
        """Divide x on y, needed for dividing freqs.

        Parameters
        ----------
        x
            The dividend
        y
            The divisor
        Returns result or x.
        """
        try:
            return float(x) / y
        except ValueError:
            return x
        except TypeError:
            return x
        except ZeroDivisionError:
            return x

    def _map_zigosity(self):
        """Map HET and HOM in ZIGOSIS columns."""
        if "HOM" in self.columns:
            self["HOM"] = self["HOM"].replace(
                {True: "HOM", np.nan: "HET", None: "HET"}
            )
            self.drop(columns="HET", inplace=True)
            self.rename(columns={"HOM": "ZIGOSITY"}, inplace=True)
        return self

    def _process_ESP6500(self):
        """Calculate frequencies on ESP6500 populations and POLYPHEN Score."""
        if "ESP6500_MAF" in self.columns:
            self[
                ["ESP6500_MAF_EA", "ESP6500_MAF_AA", "ESP6500_MAF_ALL"]
            ] = self["ESP6500_MAF"].str.split(",", expand=True)
            self["ESP6500_MAF_EA"] = self["ESP6500_MAF_EA"].apply(
                self._divide, args=(100,)
            )
            self["ESP6500_MAF_AA"] = self["ESP6500_MAF_AA"].apply(
                self._divide, args=(100,)
            )
            self["ESP6500_MAF_ALL"] = self["ESP6500_MAF_ALL"].apply(
                self._divide, args=(100,)
            )
            self.drop(columns=["ESP6500_MAF"], inplace=True)
        if "ESP6500_PH" in self.columns:
            self[["POLYPHEN_PRED", "POLYPHEN_SCORE"]] = self[
                "ESP6500_PH"
            ].str.split(":", 1, expand=True)
            self["POLYPHEN_PRED"] = (
                self["POLYPHEN_PRED"].str.strip(".").str.strip(".,")
            )
            self["POLYPHEN_SCORE"] = (
                self["POLYPHEN_SCORE"].str.split(",").str[0]
            )
            self.drop(columns=["ESP6500_PH"], inplace=True)
            self.rename(
                columns={
                    "ANNOTATION": "EFFECT",
                    "ANNOTATION_IMPACT": "IMPACT",
                    "ID": "RSID",
                },
                inplace=True,
            )
        return self

    def _process_CLINVAR(self):
        """Process and translate CLINVAR data."""
        if "CLINVAR_CLNSIG" in self.columns:
            clinvartranslation = {
                "255": "other",
                "0": "Uncertain significance",
                "1": "not provided",
                "2": "Benign",
                "3": "Likely Benign",
                "4": "Likely pathogenic",
                "5": "Pathogenic",
                "6": "drug response",
                "7": "histocompatibility",
            }
            for k, v in clinvartranslation.items():
                self["CLINVAR_CLNSIG"] = self["CLINVAR_CLNSIG"].str.replace(
                    k, v
                )
        return self

    def panel(self, panel):
        """Extract variants defined in an excel file.

        Parameters
        ----------
        panel: str
           path to the panel excel file

        Returns
        -------
        panel_df: VCFDataFrame
           Dataframe consisting of filtered variants.
        """
        panel_df = pd.DataFrame()
        gene_symbol_list = _load_panel(panel)
        for gene in gene_symbol_list:
            panel_df = panel_df.append(self.loc[self["GENE_NAME"] == gene])
        panel_df = panel_df.pipe(VCFDataFrame)
        panel_df.name = self.name
        if len(panel_df) < 1:
            raise RuntimeError("Panel Result is empty")
        return panel_df
