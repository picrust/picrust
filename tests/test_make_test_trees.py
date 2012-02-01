from unittest import main, TestCase 
from cogent import LoadTree
from cogent.parse.tree import DndParser
from ..python_scripts.make_test_trees import exclude_tip, yield_test_trees,\
  make_distance_based_exclusion_fn
"""
Tests for make_test_trees.py
"""

class TestMakeTestTrees(TestCase):
    """Tests of make_test_trees.py"""

    def setUp(self):
        self.SimpleTree = \
          DndParser("((A:0.02,B:0.01)E:0.05,(C:0.01,D:0.01)F:0.05)root;")
    
        self.SimplePolytomyTree = \
                DndParser("((A:0.02,B:0.01,B_prime:0.03)E:0.05,(C:0.01,D:0.01)F:0.05)root;")
    
    def test_exclude_tip(self):
        """exclude_tip should yield a holdout tree"""

        #Try excluding tip 'B'
        test_tree = self.SimpleTree.deepcopy()
        
        obs = exclude_tip(test_tree.getNodeMatchingName('B'),test_tree)
        obs_newick = obs.getNewick(with_distances=True)
        exp_newick = "((C:0.01,D:0.01)F:0.05,A:0.07)root;" 
        self.assertEqual(obs_newick,exp_newick)
        
        #Make sure the holdout works with a polytomy
        test_tree = self.SimplePolytomyTree.deepcopy()
        
        obs = exclude_tip(test_tree.getNodeMatchingName('B'),test_tree)
        obs_newick = obs.getNewick(with_distances=True)
        exp_newick = "((A:0.02,'B_prime':0.03)E:0.05,(C:0.01,D:0.01)F:0.05)root;" 
        self.assertEqual(obs_newick,exp_newick)

        #Make sure we raise if the tip is invalid        
         
        test_tree = self.SimpleTree.deepcopy()
        
        self.assertRaises(ValueError,exclude_tip,\
          test_tree.getNodeMatchingName('E'),test_tree)

    def test_make_distance_based_exclusion_fn(self):
        """make_distance_based_exclusion_fn should return a working function"""

        exclude_similar_strains =\
            make_distance_based_exclusion_fn(0.03)
        
        #Test that the function works
        
        test_tree = self.SimpleTree.deepcopy()
        #print test_tree.getNewick(with_distances=True)
        tip = test_tree.getNodeMatchingName('C')
        obs = exclude_similar_strains(tip,test_tree).getNewick(with_distances=True) 
        exp = "((A:0.02,B:0.01)E:0.05,C:0.06)root;"
        self.assertEqual(obs,exp)
        #Test that new function is documented
        exp_doc = 'Exclude neighbors of tip within 0.030000 branch length units'
        self.assertEqual(exp_doc,exclude_similar_strains.__doc__)


        #Test that we raise if distance is too large
        test_tree = self.SimpleTree.deepcopy()
        test_fn = make_distance_based_exclusion_fn(300.0)
        tip = test_tree.getNodeMatchingName('C')
        
        self.assertRaises(ValueError,test_fn,tip,test_tree)


    def make_distance_based_randomizer(self):
        """make_distance_based_randomizer should randomize tip labels iwthin d"""

        tip_randomizer =\
            make_distance_based_tip_label_randomizer(0.08)

        #Test that the function works

        test_tree = self.SimpleTree.deepcopy()
        print test_tree.asciiArt()
        tip = test_tree.getNodeMatchingName('C')
        obs = tip_randomizer(tip,test_tree).getNewick(with_distances=True) 
        print "OBS:",obs.asciiArt()        
        # TODO: finish this test!


    def test_yield_test_trees(self):
        """yield_test_trees should yield modified test trees"""
        
        # Test with simple tree and exclude_tip
        start_tree = self.SimpleTree.deepcopy()

        obs = yield_test_trees(start_tree, exclude_tip)
        test_trees = [tree for tree in obs]
         
        #Test that each tree  excludes the correct tip
        for i,exp_tip in enumerate(start_tree.tips()):
            node_names = [obs_tip.Name for obs_tip in test_trees[i].tips()]
            self.assertTrue(exp_tip not in node_names)
            
        #Test that the topology is correct 
        self.assertEqual(test_trees[1].getNewick(with_distances=True),\
            "((C:0.01,D:0.01)F:0.05,A:0.07)root;")

if __name__ == "__main__":
    main()
