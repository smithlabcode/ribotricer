#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools

with open('README.md') as readme_file:
    readme = readme_file.read()
with open('requirements.txt') as req_file:
    requirements = req_file.read()

setuptools.setup(
    name='RiboCop',
    version='0.2.7',
    author='Wenzheng Li',
    author_email='wenzhenl@usc.edu',
    description="Python package to detect translating ORF from Ribo-seq data",
    long_description=readme,
    long_description_content_type='text/markdown',
    url='https://github.com/wenzhenl/ribocop',
    packages=setuptools.find_packages(),
    entry_points={'console_scripts': ['RiboCop=RiboCop.cli:cli']},
    install_requires=requirements,
    classifiers=[
        'License :: OSI Approved :: GNU GPL 3',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
)
