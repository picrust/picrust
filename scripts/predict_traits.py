#!/usr/bin/env python
# File created on 15 Jul 2011
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2012, The PICRUST Project"
__credits__ = ["Jesse RR Zaneveld"]
__license__ = "GPL"
__version__ = "0.1dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"



from cogent.util.option_parsing import parse_command_line_parameters, make_option
from picrust.predict_traits import assign_traits_to_tree,\
  predict_traits_from_ancestors,make_neg_exponential_weight_fn 
from cogent import LoadTree, LoadTable
script_info = {}
script_info['brief_description'] = "Given a tree and a set of known character states (observed traits and reconstructions), output predictions for unobserved character states"
script_info['script_description'] =\
"""
[Methods description goes here]



"""
script_info['script_usage'] = [("","","")]
script_info['output_description']= "Output is a table (tab-delimited or .biom) of predicted character states"
script_info['required_options'] = [\
 # Example required option
 make_option('-i','--trait_table',type="existing_filepath",\
   help='the input trait table in tab-delimited format'),\
 make_option('-t','--tree',type="existing_filepath",help='the input Newick format tree')
]
script_info['optional_options'] = [\
 make_option('-o','--output_trait_table',type="new_filepath",default='predicted_states.tsv',help='the output filepath [default: %default]')
]
script_info['version'] = __version__



def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    # Load Tree
    tree = LoadTree(opts.tree)

    # Load Trait Table
    #TODO:  This is just a placeholder....need to specify column labels
    table = LoadTable(opts.trait_table)
    
    # Decorate tree using Trait Table
    
    # Specify the attribute where we'll store the reconstructions
    trait_label = "Reconstruction"
    tree = assign_traits_to_trees(traits,tree, trait_label=trait_label)
    
    # Perform reconstructions based using ancestral states
    
    # For now, predict all tip nodes.
    nodes_to_predict = tree.tips()
   
    #For now, use exponential weighting
    weight_fn = make_neg_exponential_weight_fn
   
    predictions = predict_traits_from_ancestors(tree,nodes_to_predict,\
     trait_label=trait_label,use_self_in_prediction = True,\
     weight_fn =weight_fn

     #TODO: Write predictions (possibly + already known traits)
     # to user specified outfile


if __name__ == "__main__":
    main()
