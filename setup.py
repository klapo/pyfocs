# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from os import path

# Get a relative path
here = path.abspath(path.dirname(__file__))


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='pyfocs',
    author="Karl Lapo and Anita Freundorfer",
    author_email='karl.lapo@uni-bayreuth.de',
    description='Processing of meteorological FODS data.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords='Fiber Optics Distributed Sensing DTS',
    url='https://github.com/klapo/pyfocs',
    packages=find_packages('src'),
    include_package_data=True,
    zip_safe=False,
    version='0.1.4.0',
    scripts=['PyFOX.py'],
    install_requires=['netcdf4',
                      'pandas',
                      'numpy',
                      'xarray<0.13,>=0.11',
                      'xmltodict',
                      'pyyaml>=5.1',
                      'dirsync',
                      'scipy',
                      'matplotlib>3'
                      ],
    python_requires='>=3.5',
    classifiers=[ "Programming Language :: Python :: 3",
                  "Development Status :: 3 - Alpha"
                  "License :: MIT License",
                  "Operating System :: OS Independent",
                  ],
    license='MIT',
    entry_points={
        'console_scripts': [
            'PyFOX.py=PyFOX.py:PyFOX.py',
        ],
    },
)
