#!/usr/bin/env python

from cogent.core.tree import TreeNode, TreeError
from cogent.parse.tree import DndParser
from cogent.util.unit_test import TestCase, main
from picrust.picrustnode import PicrustNode

class PicrustNodeTests(TestCase):
    def setUp(self):
        pass

    def test_getsubtree(self):
        """testing getting a subtree
        
        credit Julia Goodrich
        """
        otu_names = ['NineBande', 'Mouse', 'HowlerMon', 'DogFaced']
        newick = '(((Human,HowlerMon),Mouse),NineBande,DogFaced);'
        newick_reduced = '((HowlerMon,Mouse),NineBande,DogFaced);'
        tree = DndParser(newick, constructor = PicrustNode) 

        subtree = tree.unrooted().getSubTree(otu_names)
        new_tree = DndParser(newick_reduced, constructor = PicrustNode).unrooted()       
        # check we get the same names
        self.assertEqual(*[len(t.Children) for t in (subtree,new_tree)])
        self.assertEqual(subtree.getNewick(), new_tree.getNewick())


    def test_multifurcating(self):
        """Coerces nodes to have <= n children"""
        t_str = "((a:1,b:2,c:3)d:4,(e:5,f:6,g:7)h:8,(i:9,j:10,k:11)l:12)m:14;"
        t = DndParser(t_str,constructor=PicrustNode)

        # can't break up easily... sorry 80char
        exp_str = "((a:1.0,(b:2.0,c:3.0):0.0)d:4.0,((e:5.0,(f:6.0,g:7.0):0.0)h:8.0,(i:9.0,(j:10.0,k:11.0):0.0)l:12.0):0.0)m:14.0;"
        obs = t.multifurcating(2)
        self.assertEqual(obs.getNewick(with_distances=True), exp_str)
        self.assertNotEqual(t.getNewick(with_distances=True),
                            obs.getNewick(with_distances=True))

        obs = t.multifurcating(2, 0.5) 
        exp_str = "((a:1.0,(b:2.0,c:3.0):0.5)d:4.0,((e:5.0,(f:6.0,g:7.0):0.5)h:8.0,(i:9.0,(j:10.0,k:11.0):0.5)l:12.0):0.5)m:14.0;"
        self.assertEqual(obs.getNewick(with_distances=True), exp_str)

        t_str = "((a,b,c)d,(e,f,g)h,(i,j,k)l)m;"
        exp_str = "((a,(b,c))d,((e,(f,g))h,(i,(j,k))l))m;"
        t = DndParser(t_str, constructor=TreeNode)
        obs = t.multifurcating(2)
        self.assertEqual(obs.getNewick(with_distances=True), exp_str)
        obs = t.multifurcating(2, eps=10) # no effect on TreeNode type
        self.assertEqual(obs.getNewick(with_distances=True), exp_str)

        self.assertRaises(TreeError, t.multifurcating, 1)

    def test_bifurcating(self):
        """Coerces nodes to have <= 2 children"""
        t_str = "((a:1,b:2,c:3)d:4,(e:5,f:6,g:7)h:8,(i:9,j:10,k:11)l:12)m:14;"
        t = DndParser(t_str)
     
        # can't break up easily... sorry 80char
        exp_str = "((a:1.0,(b:2.0,c:3.0):0.0)d:4.0,((e:5.0,(f:6.0,g:7.0):0.0)h:8.0,(i:9.0,(j:10.0,k:11.0):0.0)l:12.0):0.0)m:14.0;"
        obs = t.bifurcating()

if __name__ == '__main__':
    main()
