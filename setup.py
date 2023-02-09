#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Part of ribotricer software
#
# Copyright (C) 2020 Saket Choudhary, Wenzheng Li, and Andrew D Smith
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import setuptools

with open("README.md") as readme_file:
    readme = readme_file.read()
with open("requirements.txt") as req_file:
    requirements = req_file.read()

setuptools.setup(
    name="ribotricer",
    version="1.3.3",
    author="Saket Choudhary, Wenzheng Li",
    author_email="saketkc@gmail.com",
    maintainer="Saket Choudhary",
    maintainer_email="saketkc@gmail.com",
    description="Python package to detect translating ORFs from Ribo-seq data",
    license="GPLv3",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/smithlabcode/ribotricer",
    packages=setuptools.find_packages(),
    entry_points={"console_scripts": ["ribotricer=ribotricer.cli:cli"]},
    python_requires=">=3.7",
    install_requires=requirements,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Utilities",
    ],
)
