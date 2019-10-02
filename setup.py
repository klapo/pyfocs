from setuptools import setup, find_packages
setup(
    name='btmm_process',
    author="Karl Lapo and Anita Freundorfer",
    author_email='karl-lapo@uni-bayreuth.de',
    description='Processing of meteorological FODS data.',
    packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
    version='0.1',
    install_requires=['netcdf4', 'pandas', 'numpy', 'xarray<0.13', 'xmltodict',
                      'pyyaml', 'dirsync', 'scipy'],
    license='MIT',
    zip_safe=False,
    include_package_data=True,
)
