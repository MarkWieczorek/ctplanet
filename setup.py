#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import versioneer

versioneer.versionfile_source = 'pycrust/_version.py'
versioneer.versionfile_build = 'pycrust/_version.py'
versioneer.tag_prefix = ''
versioneer.parentdir_prefix = 'pycrust-'

# Convert markdown README.md to restructured text (.rst) for PyPi
try:
    import pypandoc
    rst = pypandoc.convert_file('README.md', 'rst')
    long_description = rst.split('\n', 5)[5]
except(IOError, ImportError):
    print('*** pypandoc is not installed. PYPI description will not be '
          'formatted correctly. ***')
    long_description = open('README.md').read()

install_requires = ['pyshtools>=4.7']

setup(name='pycrust',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Create a crustal thickness map of a planet',
      long_description=long_description,
      url='https://github.com/MarkWieczorek/pyCrust',
      author='Mark A. Wieczorek',
      author_email='mark.a.wieczorek@gmail.com',
      license='BSD',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'Intended Audience :: Developers',
          'License :: OSI Approved :: BSD License',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Physics',
          'Topic :: Scientific/Engineering'
      ],
      keywords=['crust', 'gravity', 'geophysics'],
      packages=find_packages(),
      include_package_data=True,
      install_requires=install_requires,
      python_requires='>=3.5')
