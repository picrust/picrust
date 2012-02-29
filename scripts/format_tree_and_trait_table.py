#!/usr/bin/env python
# File created on 15 Jul 2011
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011, The PICRUST Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from os.path import splitext
from cogent import LoadTree
from cogent.util.option_parsing import parse_command_line_parameters,\
    make_option
from picrust.format_tree_and_trait_table import *


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
          make_option('-m','--tree_to_trait_mapping',default=None,type="existing_filepath",help='a two-column, tab-delimited text file mapping identifiers in the tree(column 1) to identifiers in the trait table (column 2). If supplied, the identifiers in the trait table will be converted to match the identifiers in the tree. (This mapping does not need to be supplied if the tree and trait table already use a common set of identifiers.) [default: %default]'),\
          make_option('--output_tree',default=None,type="new_filepath",help='the output tree file [default: input tree file with "_formatted" added before the extension]'),\
          make_option('--output_table_fp',default=None,type="new_filepath",help='the output tree file [default: input table filepath with "_formatted" added before the extension]'),\
          make_option('--input_table_delimiter',default='tab',type="choice",choices=delimiter_choices,\
            help='The character delimiting fields in the input trait table. Valid choices are:'+','.join(delimiter_choices)+' [default: %default]'),\
          make_option('--output_table_delimiter',default='tab',type="choice",choices=delimiter_choices,\
            help='The character delimiting fields in the output trait table. Valid choices are:'+','.join(delimiter_choices)+' [default: %default]'),\
          make_option('--suppress_bifurcating',default=False,action="store_true",help="If set, don't ensure that tree is fully bifurcating. [default: %default]"),\
          make_option('-n','--convert_to_nexus',default=False,action="store_true",help='Convert tree to NEXUS format, including a translate block mapping tip names to numbers. [default: %default]'),\
          make_option('-c','--convert_values_to_ints',default=False,action="store_true",help='Convert the values for each character state to integers. [default: %default]'),\
          make_option('--no_minimum_branch_length',default=False,action="store_true",help="If set, don't ensure all branches have at least a small but non-zero branchlength. [default: %default]"),\
          make_option('--supress_tree_filter',default=False,action="store_true",help="If set, don't filter out tree tips that aren't listed in the trait table [default: %default]"),\
          make_option('--supress_table_filter',default=False,action="store_true",help="If set, don't filter out trait table entries that aren't listed in the tree [default: %default]"),\
          make_option('-r','--add_branch_length_to_root',default=False,action="store_true",\
                      help='Add a short branch to the root node (this is required by some phylogeny programs).  The length of the branch is determined by the min_branch_length option  [default: %default]'),\
          make_option('-l','--limit_tree_to_otus_fp',type="existing_filepath",help='Will prune the reference tree to contain only those tips that are within the given OTU table')\
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
    if opts.output_table_fp:
        output_table_fp = opts.output_table_fp  
    else:
        output_table_fp = add_to_filename(trait_table_fp,"reformatted")
    
    if opts.output_tree:
        output_tree_fp = opts.output_tree  
    else:
        output_tree_fp = add_to_filename(opts.input_tree,"reformatted")

    output_reference_tree_fp = add_to_filename(output_tree_fp,"reference")
    
    delimiter_map = {"space":" ","tab":"\t","comma":","}
    input_delimiter = delimiter_map[opts.input_table_delimiter]
    output_delimiter = delimiter_map[opts.output_table_delimiter] 

    if verbose:
        print "Running with options:"
        print "\t%s:%s" %("Tree file",tree_file)
        print "\t%s:%s" %("Trait table",trait_table_fp)
        print "\t%s:%s" %("Output tree",output_tree_fp)
        print "\t%s:%s" %("Output reference tree",output_reference_tree_fp)
        print "\t%s:%s" %("Output trait table",output_table_fp)
        print "\t%s:%s" %("Add branch length to root",opts.add_branch_length_to_root)
        print "\t%s:%s" %("Convert to NEXUS?",opts.convert_to_nexus)
        print "\t%s:%s" %("Input trait table delimiter",opts.input_table_delimiter)
        print "\t%s:%s" %("Output trait table delimiter",opts.output_table_delimiter)
         
    

    # Begin reformatting
    
    
    
    root_name = "root"
    #format_for_bayestraits = True 
    #TODO: this will become a new function in the bayestraits app controller
    #if format_for_bayestraits:
    #    convert_to_nexus = True
    #    convert_to_bifurcating = True
    #    filter_table_by_tree_tips = True
    #    filter_tree_by_table_entries = True
    #    enforce_min_branch_length = True
    #    convert_trait_floats_to_ints = True
    
     
    if opts.no_minimum_branch_length:
        min_branch_length = None
    else:
        min_branch_length = 0.0001
    
    #Load inputs
    input_tree = LoadTree(tree_file)
    
    trait_table = open(trait_table_fp,"U")
    trait_table_lines = trait_table.readlines()
   
    #Get id mappings from mapping file
    if opts.tree_to_trait_mapping:
        mapping_file = open(opts.tree_to_trait_mapping,"U")
        
        trait_to_tree_mapping =\
          make_id_mapping_dict(parse_id_mapping_file(mapping_file))

    else:
        trait_to_tree_mapping = None

    #Make a clean copy of the reference tree before modifying it based on trait table
    reference_tree = input_tree.deepcopy()

    # Call reformatting function using specified parameters
    new_tree, new_trait_table_lines = \
       reformat_tree_and_trait_table(tree=input_tree,\
       trait_table_lines = trait_table_lines,\
       trait_to_tree_mapping = trait_to_tree_mapping,\
       input_trait_table_delimiter= input_delimiter,\
       output_trait_table_delimiter=output_delimiter,\
       filter_table_by_tree_tips=True,\
       convert_trait_floats_to_ints=False,\
       filter_tree_by_table_entries=True,\
       convert_to_bifurcating=True,\
       add_branch_length_to_root=False,\
       name_unnamed_nodes=True,\
       min_branch_length=min_branch_length,\
       verbose=opts.verbose) 

    # Call reformatting function using specified parameters 
    # to get reference tree
    
    new_reference_tree, not_useful_trait_table_lines =\
      reformat_tree_and_trait_table(\
      tree=reference_tree,\
      trait_table_lines = [],\
      trait_to_tree_mapping = None,\
      input_trait_table_delimiter= None,\
      output_trait_table_delimiter= None,\
      filter_table_by_tree_tips=False,\
      convert_trait_floats_to_ints=False,\
      filter_tree_by_table_entries=False,\
      convert_to_bifurcating=False,\
      add_branch_length_to_root=False,\
      name_unnamed_nodes=True,\
      min_branch_length=min_branch_length,\
      verbose=opts.verbose) 


    #Remove tips except those in trait table or in OTU table
    if opts.limit_tree_to_otus_fp:
        otu_table = open(opts.limit_tree_to_otus_fp,"U")
        otu_table_lines = otu_table.readlines()
        tips_to_keep = otu_table_lines + trait_table_lines
        tips_to_keep_in_tree = filter_table_by_presence_in_tree(new_reference_tree,tips_to_keep)
        new_reference_tree = filter_tree_tips_by_presence_in_table(new_reference_tree,\
          tips_to_keep_in_tree,verbose=opts.verbose)

        

    #Write results to files

    # Open output files
    output_trait_table_file = open(output_table_fp,"w+")
    output_tree_file  = open(output_tree_fp,"w+")
    output_reference_tree_file  = open(output_reference_tree_fp,"w+")
    
    #Output trait table file
    output_trait_table_file.write("\n".join(new_trait_table_lines))
    trait_table.close()
    output_trait_table_file.close()

    #Output tree file

    if opts.convert_to_nexus is True:
        lines = nexus_lines_from_tree(new_tree)
        output_tree_file.write("\n".join(map(str,lines)))
    else:
        output_tree_file.write(new_tree.getNewick(with_distances=True))

    output_tree_file.close() 


    #Output reference tree file
    output_reference_tree_file.write(new_reference_tree.getNewick(with_distances=True))
    output_reference_tree_file.close() 
    

if __name__ == "__main__":
    main()

   
