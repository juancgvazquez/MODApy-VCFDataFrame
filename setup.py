"""Setup File for MODApy-VCFDataFrame"""
# This file is part of the
#   MODApy-VCFDataFrame Project
#   (https://github.com/juancgvazquez/MODApy-VCFDataFrame).
# Copyright (c) 2020, Juan Carlos Vázquez
# License: MIT
#   Full Text: github.com/juancgvazquez/MODApy-VCFDataFrame/blob/master/LICENSE
import os
import pathlib

from setuptools import find_packages, setup

# EXTERNAL DATA
# Description
PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))

with open(PATH / "README.md", "r") as f:
    long_description = f.read()

# Version
with open(PATH / "VCFDataFrame" / "__init__.py") as fp:
    for line in fp.readlines():
        if line.startswith("__version__ = "):
            VERSION = line.split("=", 1)[-1].replace('"', "").strip()
            break

# Requirements
REQS = [
    "pandas",
    "xlrd",
    "joblib",
    "numpy",
    "cython",
    "cyvcf2",
    "matplotlib",
    "matplotlib-venn",
]
# SETUP #
setup(
    name="VCFDataFrame",
    version=VERSION,
    author=["Juan Carlos Vázquez"],
    author_email="juancgvazquez@gmail.com",
    description="Package to work with Variant Calling Format Files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/juancgvazquez/MODApy-VCFDataFrame/",
    license="MIT",
    install_requires=REQS,
    packages=find_packages(),
    keywords=["bioinformatics", "genomics", "vcf", "pandas"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Development Status :: 1 - Alpha",
        "Intended Audience :: Bioinformatics",
        "Intended Audience :: Genomics",
        "Topic :: Scientific/Engineering :: Bioinformatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
