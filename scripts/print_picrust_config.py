#!/usr/bin/env python

"""Print information on the PICRUSt installation.

For more details, see http://picrust.github.com
"""
from sys import platform, version as python_version, executable

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2012, The BIOM-Format Project"
__credits__ = ["Morgan Langille","Daniel McDonald", "Jose Clemente", "Greg Caporaso", 
               "Jai Ram Rideout", "Justin Kuczynski", "Andreas Wilke",
               "Tobias Paczian", "Rob Knight", "Folker Meyer", 
               "Sue Huse"]
__url__ = "http://picrust.github.com"
__license__ = "GPL"
__version__ = "0.9.1-dev"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"

try:
    from cogent import __version__ as pycogent_lib_version
except ImportError:
    pycogent_lib_version = "ERROR: Can't find the PyCogent library code - is it installed and in your $PYTHONPATH?"
    
try:
    from numpy import __version__ as numpy_lib_version
except ImportError:
    numpy_lib_version = "ERROR: Not installed - this is required! (This will also cause the BIOM library to not be importable.)"

try:
    from biom import __version__ as biom_lib_version
except ImportError:
    biom_lib_version = "ERROR: Can't find the BIOM library code (or numpy) - is it installed and in your $PYTHONPATH?"

try:
    from picrust import __version__ as picrust_lib_version
except ImportError:
    picrust_lib_version = "ERROR: Can't find the PICRUSt library code - is it installed and in your $PYTHONPATH?"

def get_script_version():
    return __version__

def print_picrust_config():
    system_info = [
     ("Platform", platform),
     ("Python/GCC version",python_version.replace('\n', ' ')),
     ("Python executable",executable)]
    
    max_len =  max([len(e[0]) for e in system_info])
    print "\nSystem information"
    print  "==================" 
    for v in system_info:
        print "%*s:\t%s" % (max_len,v[0],v[1])

    version_info = [
     ("NumPy version", numpy_lib_version),
     ("biom-format version", biom_lib_version),
     ("PyCogent version", pycogent_lib_version),
     ("PICRUSt version", picrust_lib_version),
     ("PICRUSt script version", get_script_version()),]
    
    max_len =  max([len(e[0]) for e in version_info])

    print "\nDependency versions"
    print  "===================" 
    for v in version_info:
        print "%*s:\t%s" % (max_len,v[0],v[1])
    print ""


if __name__ == '__main__':
    print_picrust_config()
