import ez_setup
ez_setup.use_setuptools()
from setuptools import setup

import bwameth

setup(name='bwameth',
      version=bwameth.__version__,
      description="align BS-Seq reads with bwa mem",
      packages=[''],
      author="Brent Pedersen",
      author_email="bpederse@gmail.com",
      license="MIT",
      install_requires=["toolshed"],
      long_description=open('README.md').read(),
      classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
      entry_points="[console_scripts]\nbwa-meth = bwameth:main",
  )
