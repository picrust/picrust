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
        
        # Just a placeholder for the moment

if __name__ == "__main__":
    main()
