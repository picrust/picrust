#!/usr/bin/env python
# File created on Jan 26 2012
from __future__ import division

__author__ = "Jesse RR Zaneveld"
__copyright__ = "Copyright 2012, The PICRUST Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Jesse RR Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"


from math import e
from cogent.util.option_parsing import parse_command_line_parameters, make_option
from numpy.ma import masked_object, array
from numpy import where, logical_not


def assign_traits_to_tree(traits, tree, trait_label="Reconstruction"):
    """Assign a dict of traits to a PyCogent tree
    
    traits -- a dict of traits, keyed by node names
    tree -- a PyCogent phylonode object
    trait_label -- a string defining the attribute in which
    traits will be recorded.  For example, if this is set to 'Reconstruction',
    the trait will be attached to node.Reconstruction
    """
    
    for node in tree.preorder():
        value_to_assign = traits.get(node.Name,None)
        #print "Assigning %s to node: %s" %(value_to_assign, node.Name)
        setattr(node,trait_label,value_to_assign)
    return tree


def linear_weight(d,max_d=1.0):
    return (max_d - d)/max_d

def make_neg_exponential_weight_fn(exp_base = 2.0):
    """Return function that exponentially weights exp_base by -d"""
    def neg_exponential_weight(d):
        return exp_base**(-1*d)
    
    return neg_exponential_weight

def equal_weight(d,constant_weight=1.0):
    return weight
    

def weighted_average_tip_prediction(tree, node_to_predict,\
  trait_label="Reconstruction", weight_fn=linear_weight):
    """Predict node traits, combining reconstructions with tip nodes
    
    tree -- a PyCogent PhyloNode tree object.   In this case,
      the tree nodes should be decorated with trait arrays stored in the
      attribute specified in trait_label.   By default this is
      node.Reconstruction.

    node_to_predict -- the Name of the node for which a reconstruction
      is desired.  Will be matched to the node.Name property on the tree.

    trait_label -- this is the node attribute where reconstructions/trait 
      values are stored
    
    weight_fn -- this is a function that takes a distance
    on the phylogenetic tree and returns a weight.  Examples include
    linear_weight (equals distance), equal_weight (a fixed value that 
    disregards distance), or neg_exponential_weight (neg. exponential
    weighting by branch length)
    """
    
    node = tree.getNodeMatchingName(node_to_predict) 
    #print dir(node)
    parent_node =  node.Parent
    sibling_nodes = node.siblings()
    
    most_rec_recon_anc =\
      get_most_recent_reconstructed_ancestor(node,trait_label)
    

    #STEP 1:  Infer traits, weight for most recently 
    #reconstructed ancestral node
    
    anc_traits = getattr(most_rec_recon_anc,trait_label,None)
    ancestor_distance =  parent_node.distance(most_rec_recon_anc)
    ancestor_weight = weight_fn(ancestor_distance)

    #STEP 2:  Infer Parent node traits
    if anc_traits is not None:
        prediction = array(anc_traits)*ancestor_weight
        total_weights = array([ancestor_weight]*len(prediction))
    else:
        prediction = None
        total_weights = None

    for child in parent_node.Children:
        child_traits = getattr(child,trait_label,None)
        if child_traits is None:
            continue
        
        distance_to_parent = parent_node.distance(child)
        weight = weight_fn(distance_to_parent)
        
        if prediction is None:
            # No ancestral states available
            prediction = array(child_traits)
        
        weights = array([weight]*len(prediction))
        
        if total_weights is None:
            total_weights = weight
            continue

        prediction += array(child_traits)*weight
        total_weights += weight

    prediction = prediction/total_weights
    
    # STEP 3: Predict target node given parent

    #NOTES: without probabilites, we're left just predicting
    # the parent





    #NOTES;  need to add a safe multiplication function
    # to handle missing values (Nones)

    #NOTE: this will need to be modified to accomodate
    # probability values

    return prediction








def predict_traits_from_ancestors(tree,nodes_to_predict,\
    trait_label="Reconstruction",use_self_in_prediction=True,\
    weight_fn=linear_weight):
    """Predict node traits given labelled ancestral states
    
    tree -- a PyCogent phylonode object, with each node decorated with the 
    attribute defined in trait label (e.g. node.Reconstruction = [0,1,1,0])

    nodes_to_predict -- a list of node names for which a trait 
    prediction should be generated
    
    trait_label -- a string defining the attribute in which the
    trait to be reconstructed is stored.  This attribute should
    contain a numpy array of trait values (which can be set to None if 
    not known)

    use_self_in_prediction -- if set to True, nodes that already
    have a trait value will be predicted using that trait value.
    If set to False, each node will ignore it's own traits when performing
    predictions, which can be useful for validation (otherwise a tree 
    would need to be generated in which known nodes have their data removed)

    

    """
    #result_tree = tree.deepcopy()
    results = {}
    n_traits = None    
    for node_label in nodes_to_predict:
        node_to_predict = tree.getNodeMatchingName(node_label)
        
        traits = getattr(node_to_predict,trait_label)
        
        # Do a little checking to make sure trait values either look 
        # like valid numpy arrays of equal length, are not specified
        # or are set to None.

        if traits is not None:
            if n_traits is None:
                #Check that trait length is consistent
                try:
                    n_traits = len(traits)
                except TypeError:
                    raise TypeError("Node trait values must be arrays!  Couldn't call len() on %s" % traits)

            if traits and len(traits) != n_traits:
                raise ValueError(\
                  "The number of traits in the array for node %s (%i) does not match other nodes (%i)" %(node_to_predict,len(traits),n_traits))
        
                
        if not use_self_in_prediction:
            # ignore knowledge about self without modifying tree
            traits = None 
        
        

        ancestral_states = \
          get_most_recent_ancestral_states(node_to_predict,trait_label)    
        
        #print "Nearest ancestral states for node:",\
        #  node_to_predict.Name,"=",ancestral_states 
        
        #Weighted average prediction from last reconstructed ancestor
        prediction =\
          weighted_average_tip_prediction(tree,node_to_predict.Name,\
          weight_fn = weight_fn)
        
        prediction = fill_unknown_traits(traits, ancestral_states)
        #print "Prediction:", prediction
        #
        results[node_label] = prediction
        
        #NOTE:  We may want to add a 'limit of accuracy' param:
        # tip nodes within 'limit of accuracy' are grouped together
        # for purposes of reconstruction
    return results





def fill_unknown_traits(traits_to_update, new_traits):
    """Returns  traits_to_update, with all None values replaced with the corresponding index in new_traits  
    
    traits_to_update -- a numpy array, with None values representing missing data
    new_traits -- a numpy array, with None values representing missing data
   
    returns an array that is traits_to_update, with all None values replaced with 
    the corresponding values in new_traits.

    This is useful when some, but not all traits of an organism are known.
    For example, PCR of specific genes may have fixed some traits, but you would
    like to predict the rest.

    Note -- in the special case where traits_to_update is the value None,
    new_traits will be returned.  This is done so that we don't have to guess
    or specify beforehand how many traits there will be when filling in tip nodes 
    """
    if traits_to_update is None:
        return array(new_traits)
    masked_traits_to_update = masked_object(traits_to_update,None)
    
    result = where(masked_traits_to_update.mask ==  True,\
      new_traits,masked_traits_to_update)  

    return array(result)

def get_most_recent_ancestral_states(node,trait_label):
    """Traverse ancestors of node until a reconstructed value is found
    
    node -- a PhyloNode object
    trait_label -- the trait attribute corresponding to 
    
    
    """
    
    ancestors = node.ancestors()
    for ancestor in node.ancestors():
        trait = getattr(ancestor,trait_label)
        if trait is not None:
            return trait
    # If we get through all ancestors, and no traits are found,
    # then there are no most recent ancestral states
    return None


def get_most_recent_reconstructed_ancestor(node,trait_label):
    """Traverse ancestors of node until a reconstructed value is found
    
    node -- a PhyloNode object
    trait_label -- the trait attribute corresponding to 
    
    
    """
    
    ancestors = node.ancestors()
    for ancestor in node.ancestors():
        trait = getattr(ancestor,trait_label)
        if trait is not None:
            return ancestor
    # If we get through all ancestors, and no traits are found,
    # then there are no most recent reconstructed ancestors
    return None


