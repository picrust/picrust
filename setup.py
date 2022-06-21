#!/usr/bin/env python

from setuptools import setup
from glob import glob


__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Greg Caporaso", "Daniel McDonald", "Jose Clemente"]
__license__ = "GPL"
__version__ = "1.1.4"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

long_description = ("PICRUSt: Phylogenetic Investigation of Communities by "
                    "Reconstruction of Unobserved States\n\n"
                    "Predictive functional profiling of microbial communities "
                    "using 16S rRNA marker gene sequences. Langille, M. "
                    "G.I.*; Zaneveld, J.*; Caporaso, J. G.; McDonald, D.; "
                    "Knights, D.; a Reyes, J.; Clemente, J. C.; Burkepile, D. "
                    "E.; Vega Thurber, R. L.; Knight, R.; Beiko, R. G.; and "
                    "Huttenhower, C. Nature Biotechnology, 1-10. 8 2013.")


setup(name='PICRUSt',
      version=__version__,
      description=('PICRUSt: Phylogenetic Investigation of Communities by '
                   'Reconstruction of Unobserved States'),
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='http://picrust.github.com',
      packages=['picrust'],
      scripts=glob('scripts/*py'),
      install_requires=['numpy == 1.22.0',
                        'cogent == 1.5.3',
                        'h5py == 2.7.1',
                        'biom-format == 2.1.7',
                        'future == 0.17.1'],
      package_data={'picrust':
                    ['data/*gz',
                     'support_files/jar/Count.jar',
                     'support_files/R/ace.R']},
      long_description=long_description)
