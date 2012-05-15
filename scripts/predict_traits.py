#!/usr/bin/env python
# File created on 15 Jul 2011
from __future__ import division
from warnings import warn

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2012, The PICRUST Project"
__credits__ = ["Jesse RR Zaneveld", "Morgan Langille"]
__license__ = "GPL"
__version__ = "0.1dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"


from warnings import warn
from math import e
from numpy import array
from cogent.util.option_parsing import parse_command_line_parameters, make_option
from cogent import LoadTree
from picrust.parse import parse_trait_table, extract_ids_from_table
from picrust.predict_traits import assign_traits_to_tree,\
  predict_traits_from_ancestors, update_trait_dict_from_file,\
  make_neg_exponential_weight_fn, biom_table_from_predictions 
from biom.table import table_factory
from cogent.util.table import Table
from picrust.format import format_biom_table

script_info = {}
script_info['brief_description'] = "Given a tree and a set of known character states (observed traits and reconstructions), output predictions for unobserved character states"
script_info['script_description'] =\
"""
This script produces predictions of unobserved traits given a phylogenetic tree and a table that summarizes which traits are present in which ancestral organisms.
In the most common usage, this script is used to predict which gene families are present in each OTU (Operational Taxonomic Unit; roughly equivalent to a bacterial 'species'), given a tree and a set of ancestral state reconstructions.

The output of the script is a trait prediction file, which summarizes the predicted traits of each organism of interest (by default, this is all of the organisms that are tips in the phylogenetic tree).

The prediction method works as follows:

    1.  For each terminal (tip) node where a prediction is to be performed, the algorithm through the reconstructed ancestral states, and finds the last node in the ancestry of our organism of interest for which a prediction is available

    2.  The trait for the organism is then predicted based on a branch-length weighted average of the ancestral node and it's close relatives. (This is necessary because technical limitations involving the handling of ambiguous characters in many Ancestral State Reconstruction programs prevent the parent node of the organism from being directly reconstructed in most cases.)

    The exact weight function to use can be specified from the commandline (see options below).

    In general, this approach causes the prediction to be a weighted average of the closest reconstructed ancestor, and the either reconstructed or directly observed trait value of the organism of interest's sibling node(s).   

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
make_option('-l','--limit_predictions_by_otu_table',type="existing_filepath",help='Specify a valid path to a legacy QIIME OTU table to perform predictions only for tips that are listed in the OTU table (regardless of abundance)'),\
make_option('-g','--limit_predictions_to_organisms',help='Limit predictions to specific, comma-separated organims ids. (Generally only useful for lists of < 10 organism ids, for example when performing leave-one-out cross-validation).'),\
make_option('-r','--reconstructed_trait_table',\
   type="existing_filepath",default=None,\
   help='the input trait table describing reconstructed traits (from ancestral_state_reconstruction.py) in tab-delimited format [default: %default]')
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
        #print value
        if value is None or len(value) == 0:
            continue
        trait_value_fields = list(value)
        new_fields.extend(trait_value_fields)
        #print new_fields
        new_line = sep.join(map(str,new_fields))
        lines.append("%s\n" % new_line)
    f.writelines(lines)
    f.close()


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    if opts.verbose:
        print "Loading tree from file:", opts.tree
    
    # Load Tree
    tree = LoadTree(opts.tree)

    table_headers =[]

    #load the asr trait table using the previous list of functions to order the arrays
    if opts.reconstructed_trait_table:
        table_headers,traits =\
                update_trait_dict_from_file(opts.reconstructed_trait_table)


    #load the trait table into a dict with organism names as keys and arrays as functions
    table_headers,genome_traits =\
            update_trait_dict_from_file(opts.observed_trait_table,table_headers)


    #Combine the trait tables overwriting the asr ones if they exist in the genome trait table.
    traits.update(genome_traits)
        
    # Specify the attribute where we'll store the reconstructions
    trait_label = "Reconstruction"
   
    if opts.verbose:
        print "Assigning traits to tree..."

    # Decorate tree using the traits
    tree = assign_traits_to_tree(traits,tree, trait_label=trait_label)
    
    if opts.verbose:
        print "Collecting list of nodes to predict..."

    #Start by predict all tip nodes.
    nodes_to_predict = [tip.Name for tip in tree.tips()]
    
    if opts.verbose:
        print "Found %i nodes to predict." % len(nodes_to_predict)

    if opts.limit_predictions_to_organisms:
        organism_id_str = opts.limit_predictions_to_organisms
        ok_organism_ids = organism_id_str.split(',')
        ok_organism_ids = [n.strip() for n in ok_organism_ids]

        if opts.verbose:
            print "Limiting predictions to user-specified ids:",\
              ok_organism_ids
        
        
        if not ok_organism_ids:
            raise RuntimeError(\
              "Found no valid ids in input: %s. Were comma-separated ids specified on the command line?"\
              % opts.limit_predictions_to_organisms)

        nodes_to_predict =\
          [n for n in nodes_to_predict if n in ok_organism_ids]

        if not nodes_to_predict:
            raise RuntimeError(\
              "Filtering by user-specified ids resulted in an empty set of nodes to predict.   Are the ids on the commmand-line and tree ids in the same format?  Example tree tip name: %s, example OTU id name: %s" %([tip.Name for tip in tree.tips()][0],ok_organism_ids[0]))
        
        if opts.verbose:
            print "After filtering organisms to predict by the ids specified on the commandline, %i nodes remain to be predicted" %(len(nodes_to_predict))
    
    if opts.limit_predictions_by_otu_table:
        if opts.verbose:
            print "Limiting predictions to ids in user-specified OTU table:",\
              opts.limit_predictions_by_otu_table
        otu_table = open(opts.limit_predictions_by_otu_table,"U")
        #Parse OTU table for ids
        
        otu_ids =\
          extract_ids_from_table(otu_table.readlines(),delimiter="\t")
        
        if not otu_ids:
            raise RuntimeError(\
              "Found no valid ids in input OTU table: %s.  Is the path correct?"\
              % opts.limit_predictions_by_otu_table)

        nodes_to_predict =\
          [n for n in nodes_to_predict if n in otu_ids]

        if not nodes_to_predict:
            raise RuntimeError(\
              "Filtering by OTU table resulted in an empty set of nodes to predict.   Are the OTU ids and tree ids in the same format?  Example tree tip name: %s, example OTU id name: %s" %([tip.Name for tip in tree.tips()][0],otu_ids[0]))
        
        if opts.verbose:
            print "After filtering by OTU table, %i nodes remain to be predicted" %(len(nodes_to_predict))
    #For now, use exponential weighting
    weight_fn = make_neg_exponential_weight_fn(e)
  
    if opts.verbose:
        print "Generating predictions..."

    # Perform predictions using reconstructed ancestral states
    predictions = predict_traits_from_ancestors(tree,nodes_to_predict,\
     trait_label=trait_label,use_self_in_prediction = True,\
     weight_fn =weight_fn,verbose=opts.verbose)

    if opts.verbose:
        print "Writing results to file: ",opts.output_trait_table

    #convert to biom format (and transpose)
    biom_predictions=biom_table_from_predictions(predictions,table_headers)

    #write biom table to file
    open(opts.output_trait_table,'w').write(\
     format_biom_table(biom_predictions))

if __name__ == "__main__":
    main()
