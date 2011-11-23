#!/usr/bin/env python
# File created on 23 Nov 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The PICRUST project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 
from os.path import exists, dirname, abspath
# Won't require PyCogent for right now, just so everyone can
# get this up and running easily. We will soon though.
#from cogent.util.unit_test import TestCase, main
from unittest import TestCase, main
from picrust.util import (get_picrust_project_dir,)

class UtilTests(TestCase):
    """ Tests of the picrust/util.py module """
    
    def setUp(self):
        """ Initialize variables: run before each test """
        pass
    
    def tearDown(self):
        """ Clean up: run after each test """
        pass
    
    def test_get_picrust_project_dir(self):
        """get_picrust_project_dir functions as expected"""
        
        # Do an explicit check on whether the file system containing
        # the current file is case insensitive.
        case_insensitive_filesystem = \
         exists(__file__.upper()) and exists(__file__.lower())
         
        actual = get_picrust_project_dir()
        # I base the expected here off the imported location of
        # picrust/util.py here, to handle cases where either the user has
        # PICRUST in their PYTHONPATH, or when they've installed it with
        # setup.py.
        # If util.py moves this test will fail -- that 
        # is what we want in this case, as the get_picrust_project_dir()
        # function would need to be modified.
        import picrust.util
        util_py_filepath = abspath(abspath(picrust.util.__file__))
        expected = dirname(dirname(util_py_filepath))
        
        if case_insensitive_filesystem:
            # make both lowercase if the file system is case insensitive
            actual = actual.lower()
            expected = expected.lower()
        self.assertEqual(actual,expected)


if __name__ == "__main__":
    main()