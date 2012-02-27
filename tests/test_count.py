#!/usr/bin/env python
# File created on 27 Feb 2012
from __future__ import division

__author__ = "Morgan Langille"
__copyright__ = "Copyright 2011, The PICRUST project"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from picrust.count import wagner_for_picrust, parse_wagner_parsimony_output
from cogent.app.util import get_tmp_filename
from cogent.util.misc import remove_files
from cogent import LoadTable
from cogent.util.table import Table

class CountTests(TestCase):
    """ """
    
    def setUp(self):
        """ """
        #create a tmp tree file
        self.in_tree1_fp = get_tmp_filename(prefix='CountTests',suffix='.nwk')
        self.in_tree1_file = open(self.in_tree1_fp,'w')
        self.in_tree1_file.write(in_tree1)
        self.in_tree1_file.close()

        #create a tmp trait file
        self.in_trait1_fp = get_tmp_filename(prefix='CountTests',suffix='.tsv')
        self.in_trait1_file=open(self.in_trait1_fp,'w')
        self.in_trait1_file.write(in_trait1)
        self.in_trait1_file.close()

        self.files_to_remove = [self.in_tree1_fp,self.in_trait1_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)
                       

    def test_wagner_for_picrust(self):
        """ test_wagner_for_picrust functions as expected with valid input
        """
        actual= wagner_for_picrust(self.in_tree1_fp,self.in_trait1_fp)
        expected=Table(['nodes','trait1','trait2'],[['11','1','3'],['12','2','3'],['10','5','2'],['14','5','3']])
        self.assertEqual(actual,expected)
   

in_tree1="""(((1:0.1,2:0.2)11:0.6,3:0.8)12:0.2,(4:0.3,D:0.4)10:0.5)14;"""

in_trait1="""tips	trait1	trait2
1	1	3
2	0	3
3	2	3
4	5	2
D	5	2"""

if __name__ == "__main__":
    main()
