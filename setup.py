from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='pyfocs',
    author="Karl Lapo and Anita Freundorfer",
    author_email='karl.lapo@uni-bayreuth.de',
    description='Processing of meteorological FODS data.',
    long_description=long_description,
    packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
    version='0.1.3.1',
    scripts=['PyFOX.py'],
    install_requires=['netcdf4',
                      'pandas',
                      'numpy',
                      'xarray<0.13',
                      'xmltodict',
                      'pyyaml>=5.1',
                      'dirsync',
                      'scipy',
                      ],
    python_requires='>=3.5',
    classifiers=[ "Programming Language :: Python :: 3",
                  "Development Status :: 3 - Alpha"
                  "License :: MIT License",
                  "Operating System :: OS Independent",
                  ],
    license='MIT',
    include_package_data=True,
)
