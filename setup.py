#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from os.path import join, split
from os import getcwd
from glob import glob

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Greg Caporaso", "Daniel McDonald", "Jose Clemente"]
__license__ = "GPL"
__version__ = "0.9.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

try:
    import numpy
except ImportError:
    raise ImportError, \
     ("biom cannot be found. Can't continue. Please install "
      "the numpy package (see www.numpy.org)")
    
try:
    import biom
except ImportError:
    raise ImportError, \
     ("biom cannot be found. Can't continue. Please install "
      "the biom package (see www.biom-format.org)")

try:
    import cogent
except ImportError:
    raise ImportError, \
     ("PyCogent cannot be found. Can't continue. Please install "
      "the PyCogent package (see www.pycogent.org).")

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
        package_data={'picrust':
                      ['data/*gz',
                       'support_files/jar/Count.jar',
                       'support_files/R/ace.R']},
        long_description=long_description)

