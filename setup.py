#!/usr/bin/env python

from setuptools import setup
from os.path import join, split
from os import getcwd
from glob import glob

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Greg Caporaso", "Daniel McDonald", "Jose Clemente"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

try:
    import numpy
except ImportError:
    raise ImportError, \
     ("numpy cannot be found. Can't continue. Please install "
      "the numpy package (see www.numpy.org)")

long_description = """PICRUSt: Phylogenetic Investigation of Communities by Reconstruction of Unobserved States
http://picrust.github.com

Manuscript in preparation. Please cite PICRUSt as http://picrust.github.com until publication.
"""

setup(name='PICRUSt',
        version=__version__,
        description='PICRUSt: Phylogenetic Investigation of Communities by Reconstruction of Unobserved States',
        author=__maintainer__,
        author_email=__email__,
        maintainer=__maintainer__,
        maintainer_email=__email__,
        url='http://picrust.github.com',
        packages=['picrust'],
        scripts=glob('scripts/*py'),
        install_requires= [
          'numpy == 1.5.1',
          'cogent == 1.5.3',
          'biom-format == 1.3.1',
          'Cython == 0.17'
        ],
        package_data={'picrust':
                      ['data/*gz',
                       'support_files/jar/Count.jar',
                       'support_files/R/ace.R']},
        long_description=long_description)
