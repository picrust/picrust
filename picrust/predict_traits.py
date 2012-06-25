#!/usr/bin/env python
# File created on Jan 26 2012
from __future__ import division

__author__ = "Jesse RR Zaneveld"
__copyright__ = "Copyright 2012, The PICRUST Project"
__credits__ = ["Jesse Zaneveld","Morgan Langille"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Jesse RR Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from collections import defaultdict
from math import e,sqrt
from random import choice
from cogent.util.option_parsing import parse_command_line_parameters, make_option
from numpy.ma import masked_object, array
from numpy import where, logical_not
from cogent.maths.stats.distribution import z_high
from cogent import LoadTable
from warnings import warn
from biom.table import table_factory,DenseOTUTable

def biom_table_from_predictions(predictions,trait_ids):
    organism_ids=predictions.keys()
    #data is in values (this transposes the matrix)
    data=map(list,zip(*predictions.values()))
    biom_table=table_factory(data,organism_ids,trait_ids, constructor=DenseOTUTable)
    return biom_table



def assign_traits_to_tree(traits, tree, trait_label="Reconstruction"):
    """Assign a dict of traits to a PyCogent tree
    
    traits -- a dict of traits, keyed by node names
    tree -- a PyCogent phylonode object
    trait_label -- a string defining the attribute in which
    traits will be recorded.  For example, if this is set to 'Reconstruction',
    the trait will be attached to node.Reconstruction
    """
    for node in tree.preorder():
        value_to_assign = traits.get(node.Name.strip(),None)
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
    

def thresholded_brownian_probability(start_state,var,d,min_val=0.0,increment=1.0,trait_prob_cutoff = 0.01):
    """Calculates the probability of a given discrete state given a continuous, brownian motion process
     
    start_state -- the starting, quantitiative value for the state (e.g. 0.36)
    var -- the variation (this is the Brownian motion parameter
      estimated during Ancestral State Reconstruction)    
    d -- the branch length along which the change will occur
    min_val -- the minimum value for the character (0 for gene copy number for example)
    increment -- amount to increase gene copy number (almost always 1.0)
    trait_prob_cutoff -- the value below which you no longer care about rare possibilities
    
    Assumes evolution by Brownian motion (directionless, suitable for neutral evolution)
    
    See Felsenstein, J. "Inferring Phylogenies" p. 430
    for further discussion of threshold models of quantitative
    traits.  The usage here is as a predictor, rather than
    for phylogenetic inference.
        
    This approach might model things like gradual loss of a gene
    by degradation better than a sudden gain by HGT, because
    intermediate values between integer copy numbers are allowed

    Algorithm:
        Variance along a branch is calculated from the supplied variance
        For each possible state under consideration from min_val --> inf 
        at increments of size 'increment', I calculate the probability of getting
        that discrete trait.

        This is calculated as the probability of having a continuous value
        in the range between that trait and the next one. So if min_val == 0 and 
        increment = 1.0 (defaults) then any continuous trait value between 0 and 
        1.0 is converted to a zero value, any value between 1.0 and 2.0 is converted to
        1, etc. 

        The probability for this interval is calculated using z_high, subtracting
        the probability of the lower bound from the probability of the upper bound
        to get the prob. of being somewhere on the interval.

        The algorithm continues until the probability falls below the 
        trait_prob_cutoff
    
    Once we have the probability of each possible change in copy number,
    we just add or subtract these from the start state to get a final value 
    
    """
    
       
    #First we calculate the probabilities of gaining or losing
    #some number of gene copies, then apply those to the actual copy numbers 
    trait_initial_value = start_state
    trait_variance =  var*d
    trait_std_dev = sqrt(trait_variance)
    mean_exp_change = 0.0 #By defn of brownian motion
    # Now we calculate the probability
    # of drawing each value of interest using the normal distribution
    
    i= min_val
    j = min_val+increment
    z_i = i/trait_std_dev
    z_j = j/trait_std_dev
    p = get_interval_z_prob(z_i,z_j)
    result = defaultdict(float)
    result[start_state+i]= p
    #print result
    while p > trait_prob_cutoff:
        i +=increment
        j +=increment
        z_i = i/trait_std_dev
        z_j = j/trait_std_dev
        p = get_interval_z_prob(i,j)
        
        if abs(p) >= trait_prob_cutoff:
            #print start_state+i,":",p
            result[start_state+i] = p
        
        if start_state - i >= min_val:
            #print start_state-i,":",p
            result[start_state - i] = p
        else:
            result[min_val] += p 
            #Values below zero should actually
            #increase the zero probability,
            #since we are using this function
            #only as an endpoint.
    
    #print result 
    
    # Now we need to get the actual
        
    return result


def get_interval_z_prob(low_z,high_z):
    """Get the probability of a range of z scores  
    """
    if low_z > high_z:
        raise ValueError(\
          "lower z value low_z must be lower than upper bound z value high_z.  Were parameters reversed?")
    #The low z represents a lower deviation from the mean,
    # and so will produce the higher probability
    high_prob = z_high(low_z)
    low_prob = z_high(high_z)
    interval_z_prob = high_prob - low_prob
    return interval_z_prob

    
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
    parent_node =  node.Parent
    #sibling_nodes = node.siblings()
    #print "Node:",node
    #print "Trait label:", trait_label 
    #print "Trait label:", trait_label 
    most_rec_recon_anc =\
      get_most_recent_reconstructed_ancestor(node,trait_label)
    
    #Preparation
    # To handle empty (None) values, we fill unknown values
    # in the ancestor with values in the tips, and vice versa


    #STEP 1:  Infer traits, weight for most recently 
    #reconstructed ancestral node
    
    if most_rec_recon_anc is not None: 
        anc_traits = getattr(most_rec_recon_anc,trait_label,None)
        ancestor_distance =  parent_node.distance(most_rec_recon_anc)
        ancestor_weight = weight_fn(ancestor_distance)
    else:
        anc_traits = ancestor_distance = ancestor_weight = None
    
    #print "Ancestor name:",most_rec_recon_anc 
    #print "Ancestor distance:",ancestor_distance 
    #print "Ancestor traits:",anc_traits

    #STEP 2:  Infer Parent node traits
    
    #TODO: abstract Step 2 to its own fn
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
    
    if prediction is not None:
        prediction = prediction/total_weights
    else:
        return None


    # STEP 3: Predict target node given parent

    #Without probabilites, we're left just predicting
    # the parent


    return prediction

def predict_random_annotated_neighbor(tree,nodes_to_predict,\
        trait_label="Reconstruction",use_self=True,verbose=False):
    """Predict traits by selecting a random, annotated tip
    tree-- PhyloNode tree object, decorated with traits (see trait_label)

    trait_label -- the attribute where arrays of reconstructed traits
    are stored.  That is, if the label is 'Reconstruction', then 
    node.Reconstruction or getattr(node,Reconstruction) should return
    an array of traits.

    use_self -- if True, allow for random prediction of self as one 
    possible outcome.

    verbose -- print verbose output.
    """

    results = {}
    n_traits = None
    annotated_nodes = \
      [t for t in tree.tips() if \
       getattr(t,trait_label,None) is not None]
    
    for node_label in nodes_to_predict:
        
        if verbose:
            print "Predicting traits for node:",node_label
        
        node_to_predict = tree.getNodeMatchingName(node_label)
        if not use_self:
            possible_nodes = [t for t in annotated_nodes if \
                t.Name != node_label]
        else:
            possible_nodes = annotated_nodes

        node_to_predict = choice(possible_nodes)
        
        if verbose:
            print "Predicting using node:", node_to_predict.Name
        
        prediction = getattr(node_to_predict,trait_label)
        results[node_label] = prediction
    return results

def predict_nearest_neighbor(tree,nodes_to_predict,\
  trait_label="Reconstruction",use_self_in_prediction=True,\
  tips_only = True, verbose = False):
    """Predict node traits given labeled ancestral states
    
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

    verbose -- output verbose debugging info 

    """
    closest_annotated_node = None
    results = {}
    n_traits = None    
    for node_label in nodes_to_predict:
        if verbose:
            print "Predicting traits for node:",node_label
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
                  "The number of traits in the array for node %s (%i) does not match other nodes (%i)" %(\
                   node_to_predict,len(traits),n_traits))
        
                
        if not use_self_in_prediction:
            # ignore knowledge about self without modifying tree
            traits = None 
       
        nearest_annotated_neighbor =\
          get_nearest_annotated_neighbor(tree,node_label,\
          trait_label=trait_label, tips_only = tips_only,\
          include_self = use_self_in_prediction)
        #print "NAN:", nearest_annotated_neighbor 
        if nearest_annotated_neighbor is None:
            raise ValueError("Couldn't find an annotated nearest neighbor for node %s on tree" % node_label)

        results[node_label] = getattr(nearest_annotated_neighbor,trait_label)
    return results

def get_nearest_annotated_neighbor(tree,node_name,\
    trait_label="Reconstruction",tips_only= True, include_self=True):
    """Return the nearest annotated node, and its distance
    
    tree -- PhyloNode object, decorated with traits in the
    attribute specified in trait_label
    node -- name of the node of interest
    trait_label -- attribute where traits are stored where 
    available
    tips_only -- if True, consider only extant, tip nodes 
    as neighbors.  if False, allow the nearest neighbor to be
    ancestral.
    """
    n1 = tree.getNodeMatchingName(node_name)
    min_dist = 99999999999999999.0
    curr_best_match = None
    if tips_only:
        neighbors = tree.tips()
    else:
        neighbors = tree.preorder()
    #print neighbors
    for n2 in neighbors:
        traits = getattr(n2,trait_label,None)
        #print n2.Name, traits
        if not traits:
            continue
        if not include_self and n1.Name == n2.Name:
            continue
        dist = n1.distance(n2) 
        if dist < min_dist:
            curr_best_match = n2
            min_dist = dist
    return curr_best_match

        
def predict_traits_from_ancestors(tree,nodes_to_predict,\
    trait_label="Reconstruction",use_self_in_prediction=True,\
    weight_fn=linear_weight, verbose = False):
    """Predict node traits given labeled ancestral states
    
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

    verbose -- output verbose debugging info 

    """
    #result_tree = tree.deepcopy()
    results = {}
    n_traits = None    
    for node_label in nodes_to_predict:
        if verbose:
            print "Predicting traits for node:",node_label
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
                  "The number of traits in the array for node %s (%i) does not match other nodes (%i)" %(\
                   node_to_predict,len(traits),n_traits))
        
                
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
         
        #prediction = fill_unknown_traits(traits, ancestral_states)
        #print "Prediction:", prediction
        #
        results[node_label] = prediction
        if verbose:
            print "First 80:",list(prediction[:min(len(prediction),80)])
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

def update_trait_dict_from_file(table_file, header = [],input_sep="\t"):
    """Update a trait dictionary from a table file

    table_file --  File name of a trait table.
    
    The first line should be a header line, with column headers equal to trait 
    (e.g. gene family) names, while the row headers should be organism 
    ids that match the tree.

    trait_dict -- a dictionary of traits, keyed by organism.  
    Items in trait dict will be overwritten if present.
    """ 
    #First line should be headers
    table=LoadTable(filename=table_file,header=True,sep=input_sep)

    #do some extra stuff to match columns if a header is provided
    if header:
        #error checking to make sure traits in ASR table are a subset of traits in genome table
        if set(header) != set(table.Header[1:]):
            if set(header).issubset(set(table.Header[1:])):
                diff_traits = set(table.Header[1:]).difference(set(header))
                warn("Missing traits in given ASR table with labels:{0}. Predictions will not be produced for these traits.".format(list(diff_traits))) 
            else:
                raise RuntimeError("Given ASR trait table contains one or more traits that do not exist in given genome trait table. Predictions can not be made.")
            
        #Note: keep the first column heading at the beginning not sorted (this is the name for the row ids
        sorted_header=[table.Header[0]]
        sorted_header.extend(header)
        table = table.getColumns(sorted_header)
    
    traits = {}
    for fields in table: 
        try:
            traits[fields[0]] = map(float,fields[1:])
        except ValueError:
            err_str =\
                    "Could not convert trait table fields:'%s' to float" %(fields[1:])
            raise ValueError(err_str)
       
    return table.Header[1:],traits
