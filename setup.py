#!/usr/bin/env python

from setuptools import setup

setup(
    name='rnapipe',
    version='0.0.1',
    author='Jessica Chung, Bernie Pope',
    author_email='jchung@unimelb.edu.au',
    packages=['rnapipe', 'rnapipe.pipeline_base'],
    entry_points={
        'console_scripts': ['rnapipe = rnapipe.main:main']
    },
    url='https://github.com/jessicachung/rnapipe',
    license='LICENSE',
    description='rnapipe is a bioinformatics pipeline for differential expression analysis\
     with support for running pipeline stages on a distributed compute cluster.',
    long_description=open('README.md').read(),
    install_requires=[
        "ruffus == 2.7.0",
        "drmaa == 0.7.6",
        "PyYAML == 3.11"
    ],
)
