#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path  # noqa: E402
from setuptools import setup, find_packages
import versioneer


install_requires = ['pyshtools>=4.7.1']

setup(name='ctplanet',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Create a crustal thickness map of a planet',
      long_description=Path('README.md').read_text(encoding='utf-8'),
      long_description_content_type='text/markdown',
      url='https://github.com/MarkWieczorek/ctplanet',
      author='Mark A. Wieczorek',
      author_email='mark.wieczorek@ipgp.fr',
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
      python_requires='>=3.6')
