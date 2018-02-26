from setuptools import setup

setup(
    name='btmm_process',
    author='Karl Lapo',
    author_email='karl.lapo@uni-bayreuth.de',
    description='Processing of meteorological data.',
    packages=['btmm_process'],
    version='0.1',
    ext_modules=cythonize('solargeo/*.pyx'),
    install_requires=['pandas', 'numpy', 'xarray'],
    tests_require=['pytest', 'engarde'])
