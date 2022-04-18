#!/usr/bin/env python

from __future__ import print_function

from setuptools import setup

try:
    from pypandoc import convert
    read_md = lambda f: convert(f, 'rst')
except ImportError:
    print('pandoc is not installed.')
    read_md = lambda f: open(f, 'r').read()

setup(name='pcalg',
      version='0.2.2',
      description='CPDAG Estimation using PC-Algorithm',
      long_description=read_md('README.md'),
      author='Keiichi SHIMA',
      author_email='keiichi@iijlab.net',
      url='https://github.com/keiichishima/pcalg/',
      py_modules=['pcalg'],
      install_requires=['networkx>=2.0', 'gsq>=0.1.6'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Information Technology',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Information Analysis',
          'Topic :: Software Development :: Libraries :: Python Modules'],
      license='BSD License',
     )
