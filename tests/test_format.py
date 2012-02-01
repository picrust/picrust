from unittest import main, TestCase 
from cogent import LoadTree
from cogent.parse.tree import DndParser
from make_test_trees import exclude_tip, yield_test_trees,\
  make_distance_based_exclusion_fn
"""
Tests for format.py
"""

class TestFormat(TestCase):
    """Tests of format.py"""

    def setUp(self):
        self.SimpleTree = \
          DndParser("((A:0.02,B:0.01)E:0.05,(C:0.01,D:0.01)F:0.05)root;")
    
        self.SimplePolytomyTree = \
                DndParser("((A:0.02,B:0.01,B_prime:0.03)E:0.05,(C:0.01,D:0.01)F:0.05)root;")
        
        # Just a placeholder for the moment

if __name__ == "__main__":
    main()
