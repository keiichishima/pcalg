#!/usr/bin/env python

from distutils.core import setup

setup(name='pcalg',
      version='0.1b1',
      description='Estimate a DAG using the pc algorithm',
      author='Keiichi SHIMA',
      author_email='keiichi@iijlab.net',
      packages=['pcalg'],
      requires=['scipy (>=0.15.1)', 'networkx (>=1.9.1)',
                'pygraphviz (>=1.2)'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Information Technology',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Information Analysis',
          'Topic :: Software Development :: Libraries :: Python Modules'],
     )
