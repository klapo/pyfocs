#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import io
import re
from glob import glob
from os import path
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup


# read the contents of your README file
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pyfocs',
    version='0.4.3',
    license='MIT',
    description='Processing of meteorological FODS data.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Karl Lapo and Anita Freundorfer',
    author_email='karl.lapo@uni-bayreuth.de',
    url='https://github.com/klapo/pyfocs',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Utilities',
    ],
    project_urls={
        'Changelog': 'https://github.com/klapo/pyfocs/blob/master/CHANGELOG.rst',
        'Issue Tracker': 'https://github.com/klapo/pyfocs/issues',
    },
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'netcdf4',
        'xarray>=0.15',
        'xmltodict',
        'pyyaml>=5.1',
        'matplotlib>3',
        'dtscalibration'
    ],

    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    # Can call the PyFOX script from the command line
    scripts=['PyFOX.py'],
)
