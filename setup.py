#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open('requirements.txt') as reqs:
    requirements = reqs.readlines()

test_requirements = requirements + ['pytest']

setup(
    name='ribocop',
    version='0.1.0',
    description="Python package to analyse ribosome profiling data",
    long_description=readme + '\n\n' + history,
    author="Saket Choudhary",
    author_email='saketkc@gmail.com',
    url='https://github.com/saketkc/ribocop',
    packages=[
        'ribocop',
    ],
    package_dir={'ribocop': 'ribocop'},
    package_data={'ribocop': ['annotation/hg38/*.*']},  # 'data/*.bw',
    #data_files=[('ribocop', ['annotation/*.*'])],
    entry_points={'console_scripts': ['ribocop=ribocop.cli:cli']},
    include_package_data=True,
    install_requires=requirements,
    license="BSD license",
    zip_safe=False,
    keywords='ribocop',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements)
