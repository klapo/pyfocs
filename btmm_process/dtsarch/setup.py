from setuptools import setup, find_packages

setup(name='dtsarch',
      version='0.0.1',
      description=u"Archiver for DTS data",
      classifiers=[],
      keywords='',
      author=u"Karl Lapo",
      author_email='karl-lapo@uni-bayreuth.de',
      url='https://github.com/klapo/dtsarch',
      license='MIT',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          'click'
      ],
      )