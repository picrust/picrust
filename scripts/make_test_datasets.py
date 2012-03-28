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

from copy import deepcopy
from cogent import LoadTree
from cogent.util.option_parsing import parse_command_line_parameters,\
  make_option
from picrust.util import make_output_dir
from picrust.parse import parse_trait_table
from picrust.format_tree_and_trait_table import filter_table_by_presence_in_tree

from picrust.make_test_datasets import yield_test_trees,\
  make_distance_based_tip_label_randomizer,make_distance_based_exclusion_fn,\
  exclude_tip,write_tree

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

    #Find which taxa are to be used in tests 
    #(by default trait table taxa)
    trait_table_header,trait_table_fields = \
            parse_trait_table(input_trait_table)


    trait_table_fields = filter_table_by_presence_in_tree(tree,\
                  trait_table_fields)


    orgs_in_table = [f[0] for f in trait_table_fields]

    #Generate a test data generation function
    #using the specified method and distance parameter
    distance_increment = 0.01
    distance_max = 0.50
    dists = [distance_increment*i for i in range(1,int(distance_max/distance_increment))]
    if opts.verbose:
        print "Building test trees for each of the following distances: %s" % (",".join(map(str,dists)))
    
    
    for curr_dist in dists:
        test_fn = method_fns[opts.method](curr_dist)



        #For now assume we want to generate tests for
        #all organisms in the table
        #TODO: add options for just doing n of these

        tips_to_predict = orgs_in_table


        base_name = "_".join(map(str,["test_tree",opts.method,curr_dist]))
        #For each tip, generate a test dataset
        test_data = []
        total_organisms = len(tips_to_predict)
        for i,tip_to_predict in enumerate(tips_to_predict):

            if opts.verbose:
                print "Generating test dataset for %i/%i: %s" %(i,total_organisms,tip_to_predict)
        
            tree_copy = tree.deepcopy()
            test_tree = \
            test_fn(tree_copy.getNodeMatchingName(tip_to_predict),tree_copy) 
        
        
            curr_filepath = write_tree(opts.output_dir,base_name,test_tree,tip_to_predict)
        
            if opts.verbose:
                print "Wrote test tree to: %s" % curr_filepath


         
    #    print "\n"
    #    print "Original Tree:\n", tree.asciiArt()
    #    print "\n"
    #    print "Test Tree:\n",test_tree.asciiArt()


if __name__ == "__main__":
    main()
