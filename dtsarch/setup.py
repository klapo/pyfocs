from setuptools import setup, find_packages
setup(
    name='btmm_process',
    author=u"Karl Lapo",
    author_email='karl-lapo@uni-bayreuth.de',
 

    description='Processing of meteorological data.',
    packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
    version='0.1',
    install_requires=['pandas', 'numpy', 'xarray'])
    license='MIT',
    zip_safe=False,
    include_package_data=True,
    install_requires=[
          'click'
          ],
)
