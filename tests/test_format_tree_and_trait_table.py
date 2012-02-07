from unittest import main, TestCase 
from cogent import LoadTree
from cogent.parse.tree import DndParser
from picrust.format_tree_and_trait_table import nexus_lines_from_tree,\
  add_branch_length_to_root,set_min_branch_length,make_nexus_trees_block,\
  filter_table_by_presence_in_tree,convert_trait_values,\
  yield_trait_table_fields,ensure_root_is_bifurcating,\
  filter_tree_tips_by_presence_in_table,print_node_summary_table,\
  add_to_filename
"""
Tests for format_tree_and_trait_tables.py
"""

class TestFormat(TestCase):
    """Tests of format.py"""

    def setUp(self):
        self.SimpleTree = \
          DndParser("((A:0.02,B:0.01)E:0.05,(C:0.01,D:0.01)F:0.05)root;")
    
        self.SimplePolytomyTree = \
                DndParser("((A:0.02,B:0.01,B_prime:0.03)E:0.05,(C:0.01,D:0.01)F:0.05)root;")
        
    def test_nexus_lines_from_tree(self):
        """Nexus lines from tree should return NEXUS formatted lines..."""
        pass

    def test_add_branch_length_to_root(self):
        """Add branch length to root should add epsilon branch lengths"""
        pass

    def test_set_min_branch_length(self):
        """set_min_branch_length should set a minimum branch length"""
        pass

    def test_make_nexus_trees_block(self):
        """make_nexus_trees_block should output a trees block in NEXUS format"""
        pass

    def test_filter_table_by_presence_in_tree(self):
        """filter_trait_table_by_presence_in_tree should filter trait table"""
        pass

    def test_convert_trait_values(self):
        """convert_trait_values should convert floats to ints in trait table"""
        pass

    def test_yield_trait_table_fields(self):
        """yield_trait_table_fields should successively yield trait table fields"""
        pass

    def test_ensure_root_is_bifurcating(self):
        """ensure_root_is_bifurcating ensures that the root of the tree is bifuracting"""
        pass

    def test_filter_tree_tips_by_presence_in_table(self):
        """filter_tree_tips_by_presence_in_table filters appropriately"""
        pass

    def test_print_node_summary_table(self):
        """print_node_summary_table prints a summary of tree nodes"""
        
        # Make sure it works for the simple tree
        tree = self.SimpleTree
        obs = [ l for l in print_node_summary_table(tree)]
        
        # Fields should be name, number of children, length, parent
        exp = ['A\t0\t0.02\tE', 'B\t0\t0.01\tE', 'E\t2\t0.05\troot',\
                'C\t0\t0.01\tF', 'D\t0\t0.01\tF', 'F\t2\t0.05\troot',\
                'root\t2\tNone\tNone']
        self.assertEqual(obs,exp) 

        # Now try it out with a tree containing a polytomy
        tree = self.SimplePolytomyTree
        obs = [ l for l in print_node_summary_table(tree)]
        
        # Fields should be name, number of children, length, parent
        exp = ['A\t0\t0.02\tE', 'B\t0\t0.01\tE','B_prime\t0\t0.03\tE',\
                'E\t3\t0.05\troot','C\t0\t0.01\tF', 'D\t0\t0.01\tF',\
                'F\t2\t0.05\troot', 'root\t2\tNone\tNone']
        self.assertEqual(obs,exp) 



    def test_add_to_filename(self):
        """add_to_filename adds text between name and file extension"""
        
        #Basic Test
        filename = "my_tree.newick"
        new_suffix = "formatted"
        delimiter = "_"

        obs = add_to_filename(filename,new_suffix,delimiter)
        exp = "my_tree_formatted.newick"
        
        self.assertEqual(obs,exp)

        #Make sure it works with multiple extensions
        
        # NOTE: this fn is based on splitext from os.path
        # therefore only the rightmost extension is considered.
        
        
        filename = "my_tree.this.is.a.newick"
        new_suffix = "formatted"
        delimiter = "."

        obs = add_to_filename(filename,new_suffix,delimiter)
        exp = "my_tree.this.is.a.formatted.newick"
        
        self.assertEqual(obs,exp)




if __name__ == "__main__":
    main()
