#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The PICRUST project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from picrust.parse import parse_marker_gene_copy_numbers

class ParseTests(TestCase):
    """ """
    
    def setUp(self):
        """ """
        self.counts_f1 = counts_f1.split('\n')
        self.counts_bad1 = counts_bad1.split('\n')
        self.counts_bad2 = counts_bad2.split('\n')
    
    def test_parse_marker_gene_copy_numbers(self):
        """ parse_marker_gene_copy_numbers functions as expected with valid input 
        """
        
        actual = parse_marker_gene_copy_numbers(self.counts_f1,
                                                metadata_identifier='cn')
        expected = {'GG_OTU_1':{'cn':2},
                    'GG_OTU_2':{'cn':4},
                    'GG_OTU_3':{'cn':1}}
        self.assertEqual(actual,expected)
    
    def test_parse_marker_gene_copy_numbers_invalid(self):
        """ parse_marker_gene_copy_numbers raises ValueError on invalid input 
        """
        self.assertRaises(ValueError,parse_marker_gene_copy_numbers,self.counts_bad1,'cn')
        self.assertRaises(ValueError,parse_marker_gene_copy_numbers,self.counts_bad2,'cn')

counts_f1 = """GG_OTU_1	2
GG_OTU_2	4
GG_OTU_3	1"""

counts_bad1 = """GG_OTU_1	hello
GG_OTU_2	4
GG_OTU_3	1"""

counts_bad2 = """GG_OTU_1	42
GG_OTU_2	4
GG_OTU_3	0"""

if __name__ == "__main__":
    main()