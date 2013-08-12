#!/usr/bin/env python

"""make_test_datasets.py


This script generates test trees for validation of approaches to character state estimation

Supported methods include:

    -- hold out analysis: systematically remove one genome at a time from the trait table
    -- collapse/randomize tree: randomize tree structure of all tips within d
    -- distance-based tree holdouts: exclude all tips within d of a taxon on the test tree
    -- jackknifing: hold out n random genomes at a time

Note that these all generate genomic test data.  Metagenomic test datasets will be generated elsewhere.


Note that although one could modify either the tree, the trait table or both,
modifying the tree is much cheaper, since it only requires (relatively) small filtered trees be saved, rather 
than complete trait tables.
"""

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.9.2-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

#from cogent import LoadTree
from sys import argv
from copy import deepcopy
import os
from random import sample, shuffle
from picrust.format_tree_and_trait_table import filter_table_by_presence_in_tree, get_sub_tree






##############################################
# Functions that modify the tree given a tip #
##############################################

def exclude_tip(tip, tree,verbose=False):
    """Return pruned tree with tip excluded
    
    tip --  a PyCogent PhyloNode object from the tree
    tree -- a PyCogent PhyloNode object (the tree)
    
    NOTE: This function exists to support distance-based exclusion of 
    neightbors.  For holdout analysis we really want to
    exclude the organism in the trait table rather than the tree.
    """
    
    # Make sure the tip is really a tip
    if not tip.isTip():
        raise ValueError("The 'holdout' fn got a non-tip node for it's 'tip' parameter.")
    
    #Old exclusion based method...but for large trees this is quite slow
    
    # Assuming that's OK, remove the tip 
    #if tip.Parent is not None:
    #    removal_ok = tip.Parent.remove(tip)
    #else:
    #    removal_ok = False
    #if not removal_ok:
    #    raise ValueError("holdout could not remove node '%s' from tree." %(tip.Name))
    #tree.prune()
    
    #Newer subtree based method
    tips_to_keep = [t.Name for t in tree.tips() if t.Name != tip.Name]
    if verbose:
        print "Generating subtree with %i tips" %len(tips_to_keep)
    subtree = get_sub_tree(tree,tips_to_keep)
    
    return subtree

def make_distance_based_exclusion_fn(exclusion_distance):
    """Returns a function that the neighbors of a tip within d.
    exclusion_distance -- the units of branchlength from which neighbors should be excluded
    
    """
    

    def exclude_tip_neighbors(tip,tree, verbose=False):
        # check that the exclusion distance won't exclude the entire tree
        max_dist = tree.maxTipTipDistance()[0]
        tips_to_delete = tip.tipsWithinDistance(exclusion_distance)
        tips_to_delete = [t.Name for t in tips_to_delete]
        if verbose:
            print "Tips to delete:",len(tips_to_delete)
        tips_to_keep = [t.Name for t in tree.tips() if t.Name not in tips_to_delete]
    
        if float(exclusion_distance) > float(max_dist) or \
          len(tips_to_delete) == len(tree.tips()):
            
            raise ValueError("specified tree would be entirely excluded because branch lengths are too short")
        
        
        #holdout_tip = tip
        holdout_tip = tip.Name
        #We're working with the *pruned* tree so we want the holdout tip
        #removed
        if holdout_tip in tips_to_keep:
            tips_to_keep.remove(holdout_tip)
            if verbose:
                print "Will remove our tip of interest from the tree: {0}".format(holdout_tip)
        if verbose:
            print "Tips to keep:", len(tips_to_keep)
        
        subtree = get_sub_tree(tree,tips_to_keep)
        del tree
        # Then it's neighbors
        #for tip in tips_to_delete:
        #    if tip == holdout_tip:
        #        #Don't exclude the tip you're testing
        #        continue
        #    tree = exclude_tip(tip,tree)
        return subtree

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
        #print "tips_to_randomize (from tipsWithinDistance):",[t.Name for t in tips_to_randomize]    
        if tip not in tips_to_randomize:
            tips_to_randomize.append(tip)
        #print "tips_to_randomize:",[t.Name for t in tips_to_randomize]    
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
    
    tree -- a PyCogent PhyloNode object for the full tree

    tree_modifier_fn -- a function that takes a tip and a tree and returns
    a new, test tree (for example a 'leave-one-out' tree)

    """
    
    
    for tip in tree.iterTips():
        # Need to pass name only to avoid
        # modifying the original tree
        tip_name = tip.Name 
        tree_copy = tree.deepcopy()
        test_tree = \
          tree_modifier_fn(tree_copy.getNodeMatchingName(tip_name),tree_copy)
        yield tip, test_tree 

def write_tree(output_dir,base_name,test_tree,tip_to_predict,delimiter='--' ):
    """Write a test trees to files with descriptive filenames"""
    file_base = delimiter.join([base_name,tip_to_predict])
    filename=os.path.join(output_dir,file_base)
    f=open(filename,"w")
    f.write(test_tree.getNewick(with_distances=True))
    f.close()
    return filename


def trait_dict_from_fields(trait_table_fields):
    """return a dict keyed by organism with values of the trait table row"""
    
    result = {}
    for fields in trait_table_fields:
        organism = fields[0]
        result[organism] = fields

    return result

def yield_genome_test_data_by_distance(tree,trait_table_fields,\
  test_fn_factory=make_distance_based_exclusion_fn, min_dist=0.0,max_dist=0.70,\
  increment=0.03,modify_tree=True,limit_to_tips=False,verbose=True):
    """Yield successive test datasets as analysis_name,test_tree,test_trait_table, test_tip tuples
    
    tree  -- a PyCogent PhyloNode object
    
    trait_table_fields -- a list of strings from a trait table
      (the first is organism name, the rest are trait values)
    
    test_fn_factory -- a factory function that returns a new test tree generating function,
      given a distance parameter.  The test tree generating function, in turn, should
      return a new test tree, given a tip of interest and the full tree. 

    min_dist -- the minimum distance at which to start generating test
      trees.

    max_dist -- the maximum distance for test trees

    increment -- the increment by which the distance parameter should be
      increased.
    
    modify_tree -- if set to False, just return the original tree without modification.

    limit_to_tips -- takes a list of organisms(tree tips) to include in the analysis.  If this is provided, instead of making test datasets 
    for all tips, only generate test datasets for tree tips with names exactly matching an entry in the supplied list of tips.
    
    For example, to limit the datasets to just E.coli K-12 MG1655 (when using the IMG / greengenes dataset)
    pass limit_to_included_tips=['NC_000913|646311926'] 
    """

    if verbose:
    #Include only organisms present in the tree
        print "Prefiltering orgs in trait table:",len([f[0] for f in trait_table_fields])
        print "Orgs with single quotes:",len([f[0] for f in trait_table_fields if "'" in f[0]]) 
        print "Tree tips with single quotes:",len([t.Name for t in tree.iterTips() if "'" in t.Name])
    
    trait_table_fields = filter_table_by_presence_in_tree(tree,\
                  trait_table_fields)

    orgs_in_table = [f[0] for f in trait_table_fields]
    if verbose:
        print "number of organisms matching between tree and table:",len(orgs_in_table)
        #print "Filtered trait table organisms:", orgs_in_table
        #print "Tree tips:",[t.Name for t in tree.iterTips()]

    if not orgs_in_table:
        raise RuntimeError("Could not match trait table organisms to tree")
    
    table_by_org = trait_dict_from_fields(trait_table_fields)

    #Generate a test data generation function
    #using the specified method and distance parameter
    
    dists = [min_dist+i*increment for i in range(int((max_dist-min_dist) // increment))]
     
    if verbose:
        print "Building test trees for each of the following distances: %s" \
        %(",".join(map(str,dists)))
    
    
    for curr_dist in dists:
        if verbose:
            print "Generating test datset for distance: %f" % curr_dist
        
        test_fn = test_fn_factory(curr_dist)


        #Select which tips will be processed

        if limit_to_tips is False:
            tips_to_predict = orgs_in_table
            if verbose:
                print "Generating test data for all tips"
        else:
            tips_to_predict = limit_to_tips
            if verbose:
                print "Generating test data for tips:",tips_to_predict
            tips_to_predict = limit_to_tips

        #For each tip of interest, generate a test dataset
        test_data = []
        total_organisms = len(tips_to_predict)
        for i,tip_to_predict in enumerate(tips_to_predict):

            if verbose:
                print "Generating test dataset (%i/%i) for: %s" %(i+1,\
                  total_organisms,tip_to_predict)
            
            if modify_tree: 
                
                if verbose:
                    print "Copying tree..."
                
                tree_copy = tree.deepcopy()
                
                if verbose:
                    print "Locating tip of interest in tree..."

                tip_of_interest = tree_copy.getNodeMatchingName(tip_to_predict)

                if verbose:
                    print "Modifying tree using function:",test_fn.__name__,\
                      test_fn.__doc__
                
                test_tree = \
                  test_fn(tip_of_interest,tree_copy) 
            else:
                test_tree = tree
            
            if verbose:
                print "Setting expected traits"
            expected_traits = table_by_org[tip_to_predict]
            
            if verbose:
                print "Generating test trait table"
            test_trait_table_fields = filter_table_by_presence_in_tree(test_tree,\
                  trait_table_fields)
            
            yield curr_dist,test_tree,tip_to_predict, expected_traits, test_trait_table_fields
        



