#!/usr/bin/env python

"""make_test_trees.py


This script generates test trees for validation of approaches to character state estimation

Supported methods include:

    -- hold out analysis: systematically remove one genome at a time from the trait table
    -- distance-based tree holdouts: exclude all tips within d of a taxon on the test tree
    -- jackknifing: hold out n random genomes at a time
    -- collapse/randomize tree: randomize tree structure of all tips within 


"""

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2007-2011, The PICRUST Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.4"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from cogent import LoadTree
from sys import argv
from copy import deepcopy
import os
from random import sample, shuffle





##############################################
# Functions that modify the tree given a tip #
##############################################

def exclude_tip(tip, tree):
    """Return pruned tree with tip excluded
    
    tip -- a PyCogent PhyloNode object from the tree
    tree -- a PyCogent PhyloNode object (the tree)
    """
    
    # Make sure the tip is really a tip
    if not tip.isTip():
        raise ValueError("The 'holdout' fn got a non-tip node for it's 'tip' parameter.")
    

    # Assuming that's OK, remove the tip 
    if tip.Parent is not None:
        removal_ok = tip.Parent.remove(tip)
    else:
        removal_ok = False
    if not removal_ok:
        raise ValueError("holdout could not remove node '%s' from tree." %(tip.Name))
    tree.prune()
    return tree

def make_distance_based_exclusion_fn(exclusion_distance):
    """Returns a function that the neighbors of a tip within d.
    exclusion_distance -- the units of branchlength from which neighbors should be excluded
    
    """
    

    def exclude_tip_neighbors(tip,tree):
        # check that the exclusion distance won't exclude the entire tree
        max_dist = tree.maxTipTipDistance()[0]
        tips_to_delete = tip.tipsWithinDistance(exclusion_distance)
        
        if float(exclusion_distance) > float(max_dist) or \
          len(tips_to_delete) == len(tree.tips()):
            
            raise ValueError("specified tree would be entirely excluded because branch lengths are too short")
        
        holdout_tip = tip
        # Then it's neighbors
        for tip in tips_to_delete:
            if tip == holdout_tip:
                #Don't exlcude the tip you're testing
                continue
            tree = exclude_tip(tip,tree)
        return tree

    new_fn = exclude_tip_neighbors
    
    #Document the new function
    new_fn_doc =\
      "Exclude neighbors of tip within %f branch length units" \
      % exclusion_distance
    new_fn.__doc__ = new_fn_doc
    
    return new_fn 


def make_distance_based_tip_label_randomizer(randomization_distance):
    """Returns a function that randomizes the labels of all nodes within d of a given node.
    randomization_distance -- the units of branchlength from which neighbors should be excluded
    """
    
    def randomize_tip_and_neighbors(tip,tree):
        
        tips_to_randomize = tip.tipsWithinDistance(randomization_distance)
        
        # Resample tips without replacement to get random labels
        new_tip_labels = [tip.Name for tip in tips_to_randomize]
        shuffle(new_tip_labels)
        
        for i,tip in enumerate(tips_to_randomize):
            tip.Name = new_tip_labels[i]
            
        return tree

    new_fn = randomize_tip_and_neighbors
    
    #Document the new function
    new_fn_doc = "Randomize tip labels within %f of tip of interest" % randomization_distance
    new_fn.__doc__ = new_fn_doc
    
    return new_fn




################################
#  Core workflow functions     #
################################


def yield_test_trees(tree,tree_modifier_fn = exclude_tip):
    """yield successive test trees from the full tree
    
    tree -- a PyCognet PhyloNode object for the full tree

    tree_modifier_fn -- a function that takes a tip and a tree and returns
    a new, test tree (for example a 'leave-one-out' tree)

    """
    
    test_set = "\t".join(["#Organism_to_predict","test_method","tree_filepath"])
    
    for tip in tree.iterTips():
        # Need to pass name only to avoid
        # modifying the original tree
        tip_name = tip.Name 
        tree_copy = tree.deepcopy()
        test_tree = \
          tree_modifier_fn(tree_copy.getNodeMatchingName(tip_name),tree_copy)
        yield tip, test_tree 

def write_trees(tips_and_test_trees, ):
    """Write test trees to files with descriptive filenames"""


def main(argv):
    """Generate test trees given parameters"""
    tree_file = argv[1]
    tree = LoadTree(tree_file)
    
    for tip_to_predict,test_tree in yield_test_trees(tree,\
      make_distance_based_tip_label_randomizer(0.15)):
        
        print "\n"
        print "Original Tree:\n", tree.asciiArt()
        print "\n"
        print "Test Tree:\n",test_tree.asciiArt() 

if __name__ == "__main__":
    main()
