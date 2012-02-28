#!/usr/bin/env python
# File created on 15 Jul 2011
from __future__ import division
from warnings import warn

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2012, The PICRUST Project"
__credits__ = ["Jesse RR Zaneveld"]
__license__ = "GPL"
__version__ = "0.1dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"


from math import e
from cogent.util.option_parsing import parse_command_line_parameters, make_option
from cogent.parse.table import SeparatorFormatParser, ConvertFields
from picrust.predict_traits import assign_traits_to_tree,\
  predict_traits_from_ancestors,make_neg_exponential_weight_fn 
from numpy import array
from cogent import LoadTree

script_info = {}
script_info['brief_description'] = "Given a tree and a set of known character states (observed traits and reconstructions), output predictions for unobserved character states"
script_info['script_description'] =\
"""
This script produces predictions of unobserved traits given a phylogenetic tree and a table that summarizes which tratis are present in which ancestral organisms.
In the most common usage, this script is used to predict which gene families are present in each OTU (Operational Taxonomic Unit; roughly equivalent to a bacterial 'species'), given a tree and a set of ancestral state reconstructions.

The output of the script is a trait prediction file, which summarizes the predicted traits of each organism of interest (by default, this is all of the organisms that are tips in the phylogenetic tree.

The prediction method works as follows:

    1.  For each terminal (tip) node where a prediction is to be performed, the algorithm through the reconstructed ancestral states, and finds the last
    node in the ancestry of our organism of interest for which a prediction is available

    2.  The trait for the organism is then predicted based on a branch-length weighted average of the ancestral node and it's close relatives.
    (This is necessary because technical limitations involving the handling of ambiguous characters in many Ancestral State Reconstruction programs
    prevent the parent node of the organism from being directly reconstructed in most cases.)

    The exact weight function to use can be specified from the commandline (see options below).

    In general, this approach causes the prediction to be a weighted average of the closest reconstructed ancestor, and the either reconstructed or directly 
    observed trait value of the organism of interest's sibling node(s).   

"""
script_info['script_usage'] = [("","","")]
script_info['output_description']= "Output is a table (tab-delimited or .biom) of predicted character states"
script_info['required_options'] = [\
make_option('-i','--observed_trait_table',type="existing_filepath",\
  help='the input trait table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format'),\
make_option('-t','--tree',type="existing_filepath",\
  help='the full reference tree, in Newick format')
]
script_info['optional_options'] = [\
 make_option('-o','--output_trait_table',type="new_filepath",\
   default='predicted_states.tsv',help='the output filepath for trait predictions [default: %default]'),\
 make_option('-r','--reconstructed_trait_table',\
   type="existing_filepath",default=None,\
   help='the input trait table describing reconstructed traits (from ancestral_state_reconstruction.py) in tab-delimited format [default: %default]')\
]
script_info['version'] = __version__

def write_results_to_file(f_out,headers,predictions,sep="\t"):
    """Write a set of predictions to a file

    headers -- a list of header columns
    predictions -- a dict of predictions, keyed by organism
    with arrays as values
    """
    f= open(f_out,"w")
    lines = [sep.join(headers)+"\n"]
    
    for pred in sorted(predictions.keys()):
        new_fields = [pred]
        value = predictions[pred]
        print value
        if value is None or len(value) == 0:
            continue
        trait_value_fields = list(value)
        new_fields.extend(trait_value_fields)
        print new_fields
        new_line = sep.join(map(str,new_fields))
        lines.append("%s\n" % new_line)
    f.writelines(lines)
    f.close()

def update_trait_dict_from_file(table_file, trait_dict = {},input_sep="\t"):
    """Update a trait dictionary from a table file

    table_file --  Lines of a trait table.
    
    The first line should be a header line, with row headers equal to trait 
    (e.g. gene family) names, while the column headers should be organism 
    ids that match the tree.

    trait_dict -- a dictionary of traits, keyed by organism.  
    Items in trait dict will be overwritten if present.
    """ 
    #First line should be headers
    table_headers = None    
    traits = trait_dict
    for i,line in enumerate(table_file):
        if i == 0:
            table_headers = line.split(input_sep)
            continue
        else:
            fields= line.split(input_sep)
            organism = fields[0]
            try:
                raw_traits = map(float,fields[1:])
            except ValueError:
                err_str =\
                  "Line %i: could not convert trait table fields:'%s' to float" %(i,fields[1:])
                raise ValueError(err_str)
            traits[organism] = raw_traits
    
    return table_headers,traits

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    # Load Tree
    tree = LoadTree(opts.tree)

    #Convert input tables to dict of traits, indexed by organism
    trait_dict ={}

    if opts.reconstructed_trait_table:
        reconstructed_trait_file = open(opts.reconstructed_trait_table,"U")
        reconstruction_table_headers,trait_dict =\
                update_trait_dict_from_file(reconstructed_trait_file)
    
    #Add in directly observed traits, overwriting ancestral state reconstructions
    #if the two overlap.
    observed_trait_file = open(opts.observed_trait_table,"U")
    observed_trait_table_headers,traits =\
            update_trait_dict_from_file(observed_trait_file,trait_dict)
    
    if observed_trait_table_headers[1:] != reconstruction_table_headers[1:]:
        warn_str =\
            "Warning!  Different headers for observed and reconstructed traits:"
        warn(warn_str)
        warn_str =\
            "Observed headers will be used in output"
        #warn(str(observed_trait_table_headers))
        #warn(str(reconstruction_table_headers))
        
    # Specify the attribute where we'll store the reconstructions
    trait_label = "Reconstruction"
   
    if opts.verbose:
        print "Assigning traits to tree..."

    # Decorate tree using the traits
    tree = assign_traits_to_tree(traits,tree, trait_label=trait_label)
    
    if opts.verbose:
        print "Collecting list of nodes to predict..."

    # For now, predict all tip nodes.
    nodes_to_predict = [tip.Name for tip in tree.tips()]
  
    #For now, use exponential weighting
    weight_fn = make_neg_exponential_weight_fn(e)
  
    if opts.verbose:
        print "Generating predictions..."

    # Perform predictions using reconstructed ancestral states
    predictions = predict_traits_from_ancestors(tree,nodes_to_predict,\
     trait_label=trait_label,use_self_in_prediction = True,\
     weight_fn =weight_fn,verbose=opts.verbose)

    if opts.verbose:
        print "Writing results to file..."
    write_results_to_file(opts.output_trait_table,observed_trait_table_headers,\
      predictions)

if __name__ == "__main__":
    main()
