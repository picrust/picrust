#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2012, The PICRUST project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"
import os
from copy import deepcopy
from cogent import LoadTree
from cogent.util.option_parsing import parse_command_line_parameters,\
  make_option
from picrust.util import make_output_dir
from picrust.parse import parse_trait_table
from picrust.format_tree_and_trait_table import filter_table_by_presence_in_tree

from picrust.make_test_datasets import yield_test_trees,\
  make_distance_based_tip_label_randomizer,make_distance_based_exclusion_fn,\
  exclude_tip,write_tree, yield_genome_test_data_by_distance

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","Generate holdout test trees from genome_tree.newick, and save results in the directory ./test_holdout_trees/.","%prog -t genome_tree.newick -o ./test_holdout_trees")]
script_info['output_description']= ""
method_choices = ['exclude_tips_by_distance','randomize_tip_labels_by_distance','collapse_tree_by_distance']
script_info['required_options'] = [
 
 make_option('-i','--input_trait_table',type='existing_filepath',\
   help='the input trait table.'),\
 make_option('-t','--input_tree',type='existing_filepath',\
   help='the input tree in Newick format'),\
]
script_info['optional_options'] = [\
  make_option('-o','--output_dir',default='./test_datasets/',type='new_dirpath',\
  help='the output directory.  Duplicate trees, trait tables, expected values and prediction files will be saved here.[default:%default]'),\
  make_option('-m','--method',type='choice',\
    choices=method_choices,default=method_choices[0],\
    help='The test method to use in generating test data.  Valid choices are:'\
      +','.join(method_choices)+' [default: %default]'),\
  
]

script_info['version'] = __version__


    
def main():
    """Generate test trees given parameters"""
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    

    input_trait_table = open(opts.input_trait_table,"U")
    tree = LoadTree(opts.input_tree)
    make_output_dir(opts.output_dir)
    

    method_fns =\
      {"exclude_tips_by_distance":\
         make_distance_based_exclusion_fn,\
       "randomize_tip_labels_by_distance":\
         make_distance_based_tip_label_randomizer
       }

    test_fn_factory = method_fns[opts.method]
    
    #Find which taxa are to be used in tests 
    #(by default trait table taxa)
    trait_table_header,trait_table_fields = \
            parse_trait_table(input_trait_table)


    trait_table_fields = [t for t in trait_table_fields]
    
    test_datasets = \
      yield_genome_test_data_by_distance(tree,trait_table_fields,test_fn_factory,\
      min_dist = 0.0, max_dist=0.45,increment=0.03, verbose = opts.verbose)
    
    
    for curr_dist,test_tree,tip_to_predict,expected_traits in test_datasets:    
       

        #Write tree
        base_name = "--".join(map(str,["test_tree",opts.method,curr_dist]))
        curr_filepath = write_tree(opts.output_dir,base_name,test_tree,tip_to_predict)
        if opts.verbose:
            print "Wrote test tree to: %s" % curr_filepath
        
        #Write expected trait table
        base_name = "--".join(map(str,["exp_traits",opts.method,curr_dist,tip_to_predict]))
                
        exp_trait_table_lines = [trait_table_header]
        exp_trait_table_lines.append("\t".join(expected_traits)+"\n")
        print "Expected_trait_table_lines:",exp_trait_table_lines
        filename=os.path.join(opts.output_dir,base_name)
        
        if opts.verbose:
            print "Writing expected trait table to:", filename
        
        f=open(filename,"w")
        f.write("".join(exp_trait_table_lines))
        f.close()

        #Write test trait table
        test_trait_table_fields = deepcopy(trait_table_fields)
        test_trait_table_fields.remove(expected_traits)
        test_trait_table_lines = [trait_table_header]
        test_trait_table_lines.extend(["\t".join(r)+"\n" for r in test_trait_table_fields])
        
        #print "Test_trait_table_lines:",test_trait_table_lines
        base_name = "--".join(map(str,["test_trait_table",opts.method,curr_dist,tip_to_predict]))
        filename=os.path.join(opts.output_dir,base_name)
        
        if opts.verbose:
            print "Writing test trait table to:", filename
        
        f=open(filename,"w")
        f.write("".join(test_trait_table_lines))
        f.close()




         
    #    print "\n"
    #    print "Original Tree:\n", tree.asciiArt()
    #    print "\n"
    #    print "Test Tree:\n",test_tree.asciiArt()


if __name__ == "__main__":
    main()
