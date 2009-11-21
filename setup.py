#!/usr/bin/python

from distutils.core import setup

classifiers=[
      'Development Status :: 4 - Beta',
      'Environment :: X11 Applications',
      'Intended Audience :: End Users/Desktop',
      'License :: OSI Approved :: GNU General Public License (GPL)',
      'Operating System :: POSIX :: Linux',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering',
    ]


setup(name='pyraman',
      version='0.2.1',
      description = 'Lib to process Raman spectra in pylab',
      author = 'Alexey Brazhe',
      py_modules = ['raman'],
      classifiers=classifiers,
      )

