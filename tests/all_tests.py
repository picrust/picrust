#!/usr/bin/env python
#
# suite of biom format unit tests.
# run suite by executing this file
# addapted from PyCogent and biom-format alltests.py
# https://github.com/pycogent
# https://github.com/biom-format/

import sys, os
from os.path import split, splitext, join
from glob import glob
import cogent.util.unit_test as unittest

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Jose Carlos Clemente Litran", "Daniel McDonald",
               "Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

def my_import(name):
    """Imports a module, possibly qualified with periods. Returns the module.
    
    __import__ only imports the top-level module.
    
    Recipe from python documentation at:
    http://www.python.org/doc/2.4/lib/built-in-funcs.html
    """
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod

def suite():
    test_dir = split(__file__)[0]
    modules_to_test = [splitext(split(e)[1])[0] for e in glob(join(test_dir,"test_*.py"))]

    alltests = unittest.TestSuite()
    
    for module in modules_to_test:
        test = unittest.findTestCases(my_import(module))
        alltests.addTest(test)
    return alltests

class BoobyTrappedStream(object):
    def __init__(self, output):
        self.output = output

    def write(self, text):
        self.output.write(text)
        raise RuntimeError, "Output not allowed in tests"
        
    def flush(self):
        pass
        
    def isatty(self):
        return False

if __name__ == '__main__':
    if '--debug' in sys.argv:
        s = suite()
        s.debug()
    else:
        orig = sys.stdout
        if '--output-ok' in sys.argv:
            sys.argv.remove('--output-ok')
        else:
            sys.stdout = BoobyTrappedStream(orig)
        try:
            unittest.main(defaultTest='suite', argv=sys.argv)
        finally:
            sys.stdout = orig
