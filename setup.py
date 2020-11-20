"""Setup File for MODApy-VCFDataFrame"""
# This file is part of the
#   MODApy-VCFDataFrame Project
#   (https://github.com/juancgvazquez/MODApy-VCFDataFrame).
# Copyright (c) 2020, Juan Carlos Vázquez
# License: MIT
#   Full Text: github.com/juancgvazquez/MODApy-VCFDataFrame/blob/master/LICENSE

from setuptools import find_packages, setup

# EXTERNAL DATA
# Description
with open("./README.md", "r") as f:
    long_description = f.read()
# Version
version = {}
with open("./VCFDataFrame/__init__.py", "r") as v:
    exec(v.read(), version)
# Requirements
requirements = ["cyvcf2", "numpy", "pandas"]


# SETUP #
setup(
    name="VCFDataFrame",
    version=version["__version__"],
    author=["Juan Carlos Vázquez"],
    author_email="juancgvazquez@gmail.com",
    description="Package to work with Variant Calling Format Files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/juancgvazquez/MODApy-VCFDataFrame/",
    license="MIT",
    install_requires=requirements,
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
