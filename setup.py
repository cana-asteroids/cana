#!/usr/bin/env python

from setuptools import setup

# PACKAGE METADATA
# ##################
NAME = 'cana-asteroids'
FULLNAME = "CANA"
VERSION = '0.11'
DESCRIPTION = 'Codes for ANalysis of Asteroids'
with open("README.rst") as f:
    LONG_DESCRIPTION = ''.join(f.readlines())
AUTHOR = '''M. De Pra, J. Carvano, D. Morate,
            J. Licandro, N. Pinilla-Alonso'''
AUTHOR_EMAIL = 'mariondepra@gmail.com'
MAINTAINER = 'M. De Pra'
MAINTAINER_EMAIL = AUTHOR_EMAIL
URL = 'https://github.com/depra/cana'
LICENSE = 'MIT License'

# TO BE INSTALLED
# ##################
PACKAGES = ['cana',
            'cana.spectools',
            'cana.pipelines',
            'cana.datasets',
            'cana.composition',
            'cana.thermal']

PACKAGE_DATA = {
   'cana.datasets': ['data/taxonomy/*',
                     'data/photometry/*/*',
                     'data/testdata/*/*/*']
}

# DEPENDENCIES
# ##################
INSTALL_REQUIRES = [
    'numpy',
    'pandas',
    'scipy',
    'sklearn',
    'matplotlib',
    'dataclasses',
    'numba',
    'dask[complete]',
    # 'http://github.com/PolyChord/PolyChordLite@master'
]

# DEPENDENCY_LINKS = ['git+https://github.com/PolyChord/PolyChordLite@master']

PYTHON_REQUIRES = ">=3.6"

if __name__ == '__main__':
    setup(name=NAME,
          description=DESCRIPTION,
          #   long_description=LONG_DESCRIPTION,
          version=VERSION,
          author=AUTHOR,
          author_email=AUTHOR_EMAIL,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          license=LICENSE,
          url=URL,
          #   platforms=PLATFORMS,
          #   scripts=SCRIPTS,
          packages=PACKAGES,
           # dependency_links = DEPENDENCY_LINKS,
          #   ext_modules=EXT_MODULES,
          package_data=PACKAGE_DATA,
          #   classifiers=CLASSIFIERS,
          #   keywords=KEYWORDS,
          #   cmdclass=CMDCLASS,
          install_requires=INSTALL_REQUIRES,
          python_requires=PYTHON_REQUIRES,
          )
