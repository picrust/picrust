#!/usr/bin/env python
# File created on 15 Jul 2011
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011, The PICRUST Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from os.path import splitext
from cogent import LoadTree
from cogent.util.option_parsing import parse_command_line_parameters,\
    make_option


# Set up commandline parameters
script_info = {}
script_info['brief_description'] = "Formatting script for filtering and reformatting trees and trait tables."
script_info['script_description'] =\
  """Reformats scripts and trait tables.  Optional fixes include:  
        -- Add short (epsilon) branch lengths in place of 0 length branches 
        -- Filter out taxa that don't match between tree and trait table 
        -- Output tree in NEXUS format 
        -- Ensure tree is bifurcating (remove polytomies using very short branches)
        -- Convert floating point trait values to integers
        -- Add a short branch length to the root branch (required by BayesTraits)
        -- Remove internal node names (required by BayesTraits)
        """

script_info['script_usage'] = [("Reformat a tree and trait table with default options:",\
    "Matches taxa in tree and trait table, reformats, and saves as NEXUS",\
    "%prog -t ./example_tree.tree -i ./example_trait_table.txt")]
script_info['output_description']= "Outputs a reformatted tree and trait table."
script_info['required_options'] = [\
          make_option('-t','--input_tree',type="existing_filepath",help='the input tree (Newick format)'),\
          make_option('-i','--input_trait_table',type="existing_filepath",help='the input trait table (QIIME OTU table format)')
                  ]

delimiter_choices = ['tab','space','comma']
script_info['optional_options'] = [\
          make_option('--output_tree',default=None,type="new_filepath",help='the output tree file [default: input tree file with "_formatted" added before the extension]'),\
          make_option('--output_table',default=None,type="new_filepath",help='the output table file [default: input tree file with "_formatted" added before the extension]'),\
          make_option('--input_table_delimiter',default='tab',type="choice",choices=delimiter_choices,\
            help='The character delimiting fields in the input trait table. Valid choices are:'+','.join(delimiter_choices)+' [default: %default]'),\
          make_option('--output_table_delimiter',default='tab',type="choice",choices=delimiter_choices,\
            help='The character delimiting fields in the output trait table. Valid choices are:'+','.join(delimiter_choices)+' [default: %default]'),\
          make_option('--suppress_bifurcating',default=False,action="store_true",help="If set, don't ensure that tree is fully bifurcating. [default: %default]"),\
          make_option('-n','--convert_to_nexus',default=False,action="store_true",help='Convert tree to NEXUS format, including a translate block mapping tip names to numbers. [default: %default]'),\
          make_option('-c','--convert_values_to_ints',default=False,action="store_true",help='Convert the values for each character state to integers. [default: %default]'),\
          make_option('-m','--no_minimum_branch_length',default=False,action="store_true",help="If set, don't ensure all branches have at least a small but non-zero branchlength. [default: %default]"),\
          make_option('--supress_tree_filter',default=False,action="store_true",help="If set, don't filter out tree tips that aren't listed in the trait table [default: %default]"),\
          make_option('--supress_table_filter',default=False,action="store_true",help="If set, don't filter out trait table entries that aren't listed in the tree [default: %default]"),\
          make_option('-r','--add_branch_length_to_root',default=False,action="store_true",\
            help='Add a short branch to the root node (this is required by some phylogeny programs).  The length of the branch is determined by the min_branch_length option  [default: %default]')\
                  ]
script_info['version'] = __version__

def main():

    # Parse input to get parameters 
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)
    
    tree_file = opts.input_tree
    trait_table_fp = opts.input_trait_table
    verbose = opts.verbose 
    
    #Handle parameters with more complex defaults
    if opts.output_table:
        output_table_fp = opts.output_table  
    else:
        output_table_fp = add_to_filename(trait_table_fp,"reformatted")
    
    if opts.output_tree:
        output_tree_fp = opts.output_tree  
    else:
        output_tree_fp = add_to_filename(trait_table_fp,"reformatted")

     

    if verbose:
        print "Running with options:"
        print "%s:%s" %("Tree file",tree_file)
        print "%s:%s" %("Trait table",trait_table_fp)
        print "%s:%s" %("Output tree",output_tree_fp)
        print "%s:%s" %("Output trait table",output_table_fp)
        print "%s:%s" %("Add branch length to root",opts.add_branch_length_to_root)
        print "%s:%s" %("Convert to NEXUS?",opts.convert_to_nexus)
         
    # Open output files
    output_trait_table_file = open(output_table_fp,"w+")
    output_tree_file  = open(output_tree_fp,"w+")
    

    # Begin reformatting
    
    # TODO: convert parameter candidates to params
    #delimiter = "\t"
    
    
    delimiter = " "
    min_branch_length = 0.0001
    root_name = "root"
    
    if format_for_bayestraits:
        convert_to_nexus = True
        convert_to_bifurcating = True
        filter_table_by_tree_tips = True
        filter_tree_by_table_entries = True
        enforce_min_branch_length = True
        convert_trait_floats_to_ints = True

    #Load inputs
    input_tree = LoadTree(tree_file)
    
    trait_table = open(trait_table_fp,"U")
    trait_table_lines = trait_table.readlines()
    
    #if opts.verbose:
        #print "Tree prior to reformatting:"
        #print input_tree.getNewick(with_distances=True)
        #print "Tree nodes prior to reformatting:"
        #print_node_summary_table(input_tree)
    
    # Remove taxa missing in tree tips
    if filter_table_by_tree_tips: 
        trait_table_lines = filter_table_by_presence_in_tree(input_tree,trait_table_lines,delimiter=delimiter)
    
    #Convert floating point values to ints
    if convert_trait_floats_to_ints:
        trait_table_lines = convert_trait_values(\
            trait_table_lines,conversion_fn = int,delimiter=delimiter)
   
    #Write out results
    output_trait_table_file.writelines(trait_table_lines)
    trait_table.close()
    output_trait_table_file.close()
   
    # Use the reformatted trait table to filter tree
    trait_table = open(output_table_fp,"U")
    trait_table_lines = trait_table.readlines()
    

    if filter_tree_by_table_entries:
        input_tree = filter_tree_tips_by_presence_in_table(input_tree,trait_table_lines,delimiter=delimiter) 
    
    # Tree reformatting
    
    
    if convert_to_bifurcating:
        input_tree = input_tree.bifurcating() # Required by most ancSR programs
    
   
    if convert_to_bifurcating:
        input_tree = ensure_root_is_bifurcating(input_tree)
        # The below nutty-looking re-filtering step is necessary
        # When ensuring the root is bifurcating, internal nodes can get moved to the tips
        # So without additional filtering we get unannotated tip nodes
        if filter_tree_by_table_entries:
            input_tree = filter_tree_tips_by_presence_in_table(input_tree,trait_table_lines,delimiter=delimiter) 
    
    if enforce_min_branch_length:
        input_tree = set_min_branch_length(input_tree,min_length = min_branch_length)
    
    if opts.add_branch_length_to_root:
        input_tree = add_branch_length_to_root(input_tree,root_name=input_tree.Name,root_length=min_branch_length)
    
    input_tree.prune() 
    
    #if opts.verbose:
    #    print "Tree nodes after reformatting:"
    #    print_node_summary_table(input_tree)
    #    print "Convert to nexus is True?",convert_to_nexus is True
    if convert_to_nexus is True:
        lines = nexus_lines_from_tree(input_tree)
        output_tree_file.write("\n".join(map(str,lines)))
    else:
        output_tree_file.write(input_tree.getNewick(with_distances=True))
    
    output_tree_file.close()

if __name__ == "__main__":
    main()

   
