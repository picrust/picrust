#!/usr/bin/env python
# File created on 23 Nov 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Greg Caporaso","Morgan Langille"]
__license__ = "GPL"
__version__ = "0.9.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 
from os.path import exists, dirname, abspath
# Won't require PyCogent for right now, just so everyone can
# get this up and running easily. We will soon though.
#from cogent.util.unit_test import TestCase, main
from unittest import TestCase, main
from picrust.util import (get_picrust_project_dir,)

from cogent.core.tree import TreeError
from cogent.parse.tree import DndParser
from picrust.util import PicrustNode,\
  transpose_trait_table_fields

class PicrustNodeTests(TestCase):
    def setUp(self):
        pass

    def test_getsubtree(self):
        """testing getting a subtree
        """
        otu_names = ['NineBande', 'Mouse', 'HowlerMon', 'DogFaced']
        newick = '(((Human,HowlerMon),Mouse),NineBande,DogFaced);'
        newick_reduced = '((Mouse,HowlerMon),NineBande,DogFaced);'
        tree = DndParser(newick, constructor = PicrustNode) 

        subtree = tree.getSubTree(otu_names)
        new_tree = DndParser(newick_reduced, constructor = PicrustNode)       
        # check we get the same names
        self.assertEqual(*[len(t.Children) for t in (subtree,new_tree)])
        self.assertEqual(subtree.getNewick(), new_tree.getNewick())


#    def test_multifurcating(self):
#        """Coerces nodes to have <= n children"""
#        t_str = "((a:1,b:2,c:3)d:4,(e:5,f:6,g:7)h:8,(i:9,j:10,k:11)l:12)m:14;"
#        t = DndParser(t_str)

        # can't break up easily... sorry 80char
#        exp_str = "((a:1.0,(b:2.0,c:3.0):0.0)d:4.0,((e:5.0,(f:6.0,g:7.0):0.0)h:8.0,(i:9.0,(j:10.0,k:11.0):0.0)l:12.0):0.0)m:14.0;"
#        obs = t.multifurcating(2)
#        self.assertEqual(obs.getNewick(with_distances=True), exp_str)
#        self.assertNotEqual(t.getNewick(with_distances=True),
#                            obs.getNewick(with_distances=True))

#        obs = t.multifurcating(2, 0.5) 
#        exp_str = "((a:1.0,(b:2.0,c:3.0):0.5)d:4.0,((e:5.0,(f:6.0,g:7.0):0.5)h:8.0,(i:9.0,(j:10.0,k:11.0):0.5)l:12.0):0.5)m:14.0;"
#        self.assertEqual(obs.getNewick(with_distances=True), exp_str)

#        t_str = "((a,b,c)d,(e,f,g)h,(i,j,k)l)m;"
#        exp_str = "((a,(b,c))d,((e,(f,g))h,(i,(j,k))l))m;"
#        t = DndParser(t_str, constructor=PicrustNode)
#        obs = t.multifurcating(2)
#        self.assertEqual(obs.getNewick(with_distances=True), exp_str)
#        obs = t.multifurcating(2, eps=10) # no effect on TreeNode type
#        self.assertEqual(obs.getNewick(with_distances=True), exp_str)

#        self.assertRaises(TreeError, t.multifurcating, 1)

    def test_bifurcating(self):
        """Coerces nodes to have <= 2 children"""
        t_str = "((a:1,b:2,c:3)d:4,(e:5,f:6,g:7)h:8,(i:9,j:10,k:11)l:12)m:14;"
        t = DndParser(t_str)
     
        # can't break up easily... sorry 80char
        exp_str = "((a:1.0,(b:2.0,c:3.0):0.0)d:4.0,((e:5.0,(f:6.0,g:7.0):0.0)h:8.0,(i:9.0,(j:10.0,k:11.0):0.0)l:12.0):0.0)m:14.0;"
        obs = t.bifurcating()



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

    def test_transpose_trait_table_fields(self):
        """transpose_table_fields should correctly transpose table fields"""
        # The purpose is to make it easy to transpose tab-delimited tables 
        # before analysis starts (and to have an internal method where needed
        
        # Get header
        test_header = "genomes\t123456\t123457\t123458\t123459\n"
        test_data_fields =[\
          ["gene1",0,3,5,0],\
          ["gene2",1,2,1,0],\
          ["gene3",0,0,0,1]]

        exp_header  = "genomes\tgene1\tgene2\tgene3\n"
         
        exp_data_fields =[\
          ['123456',0,1,0],\
          ['123457',3,2,0],\
          ['123458',5,1,0],\
          ['123459',0,0,1]]
        
        new_header,new_data_fields = transpose_trait_table_fields(test_data_fields,\
          header=test_header,id_row_idx=0)
        
        self.assertEqual(new_header,exp_header)
        self.assertEqual(new_data_fields,exp_data_fields)

        
        pass


if __name__ == "__main__":
    main()
