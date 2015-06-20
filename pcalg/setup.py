#!/usr/bin/env python

from distutils.core import setup

setup(name='pcalg',
      version='0.1b1',
      description='Estimate a DAG using the pc algorithm',
      author='Keiichi SHIMA',
      author_email='keiichi@iijlab.net',
      py_modules=['pcalg'],
      requires=['networkx (>=1.9.1)', 'pygraphviz (>=1.2)',
                'gsq (>=0.1b1)'],
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
