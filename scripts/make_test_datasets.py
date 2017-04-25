#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.1.1-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"
import os
from copy import deepcopy
from cogent.parse.tree import DndParser
from cogent.util.option_parsing import parse_command_line_parameters,\
  make_option
from biom.table import Table
from picrust.util import parse_table_to_biom, make_output_dir,\
     PicrustNode, write_biom_table,transpose_trait_table_fields
from picrust.parse import parse_trait_table
from picrust.format_tree_and_trait_table import filter_table_by_presence_in_tree,\
  set_label_conversion_fns, fix_tree_labels, convert_trait_table_entries
from numpy import array
from picrust.make_test_datasets import yield_test_trees,\
  make_distance_based_tip_label_randomizer,make_distance_based_exclusion_fn,\
  exclude_tip,write_tree, yield_genome_test_data_by_distance


script_info = {}
script_info['brief_description'] = "Generates test datasets for cross-validation studies of PICRUSt's accuracy"
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
  make_option('--min_dist',default=0.0,type='float',\
  help='the minimum phylogenetic distance to use with the holdout method, if applicable.  Usually 0.0.[default:%default]'),\
   make_option('--suppress_tree_modification',default=False,action="store_true",help='If passed, modify only the trait table, not the tree . [default: %default]'),\
  make_option('--dist_increment',default=0.03,type='float',\
  help='the phylogenetic distance increment to use with the holdout method, if applicable.[default:%default]'),\
  make_option('--max_dist',default=0.45,type='float',\
  help='the maximum phylogenetic distance to use with the holdout method, if applicable.[default:%default]'),\
  make_option('--limit_to_tips',default='',type='string',\
  help='if specified, limit test dataset generation to specified tips (comma-separated).[default:%default]'),\
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

    if opts.verbose:
        print "Loading trait table..."
    input_trait_table = open(opts.input_trait_table,"U")

    if opts.verbose:
        print "Loading tree..."
    #PicrustNode seems to run into very slow/memory intentsive perfromance...
    #tree = DndParser(open(opts.input_tree),constructor=PicrustNode)
    tree = DndParser(open(opts.input_tree))

    if opts.verbose:
        print "Parsing trait table..."
    #Find which taxa are to be used in tests
    #(by default trait table taxa)
    trait_table_header,trait_table_fields = \
            parse_trait_table(input_trait_table)

    if opts.verbose:
       print "Ensuring tree and trait table labels are formatted consistently..."

    label_conversion_fns = set_label_conversion_fns(verbose=opts.verbose)

    fix_tree_labels(tree,label_conversion_fns)

    trait_table_fields = convert_trait_table_entries(trait_table_fields,\
      value_conversion_fns = [],\
      label_conversion_fns = label_conversion_fns)

    trait_table_fields = [t for t in trait_table_fields]
    print "Number of trait table fields with single quotes:",\
     len([t for t in trait_table_fields if "'" in t[0]])

    if opts.verbose:
        print "Making output directory..."
    make_output_dir(opts.output_dir)


    if opts.limit_to_tips:

        included_tips = opts.limit_to_tips.split(",")
        if opts.verbose:
            print "Limiting test datasets to %i tips: %s" %(len(included_tips),included_tips)
    else:
        included_tips = False

    method_fns =\
      {"exclude_tips_by_distance":\
         make_distance_based_exclusion_fn,\
       "randomize_tip_labels_by_distance":\
         make_distance_based_tip_label_randomizer
       }

    test_fn_factory = method_fns[opts.method]

    if opts.verbose:
        print "Setting tree modification method to:", opts.method
        print "(%s)" % test_fn_factory.__doc__

    modify_tree = True
    if opts.suppress_tree_modification:
        if opts.verbose:
            print "Suppressing modification of tree when making test datasets"
        modify_tree = False

    if opts.verbose:
        print "Starting generation of test datsets"

    test_datasets = \
      yield_genome_test_data_by_distance(tree,trait_table_fields,\
      test_fn_factory,min_dist = opts.min_dist,\
      max_dist=opts.max_dist,increment=opts.dist_increment,\
      modify_tree=modify_tree,limit_to_tips= included_tips,verbose = opts.verbose)

    if opts.verbose:
        print "Writing files for test  datasets"

    for curr_dist,test_tree,tip_to_predict,\
        expected_traits,test_trait_table_fields in test_datasets:

        if included_tips is not False:
            if tip_to_predict not in included_tips:
                if opts.verbose:
                    print "Skipping tip %s: limiting to tip(s): %s" %(tip_to_predict,included_tips)
                continue


        #Make a safe version of tip to predict
        # So odd characters like | don't mess up OS

        safe_tip_to_predict = "'%s'"%tip_to_predict

        #Write tree
        base_name = "--".join(map(str,["test_tree",opts.method,curr_dist]))
        curr_filepath = write_tree(opts.output_dir,base_name,test_tree,safe_tip_to_predict)
        if opts.verbose:
            print "Wrote test tree to: %s" % curr_filepath

        #Write expected trait table
        base_name = "--".join(map(str,["exp_traits",opts.method,curr_dist,safe_tip_to_predict]))

        exp_trait_table_lines = [trait_table_header]
        exp_trait_table_lines.append("\t".join(expected_traits)+"\n")
        #print "Expected_trait_table_lines:",exp_trait_table_lines
        filename=os.path.join(opts.output_dir,base_name)
        if opts.verbose:
            print "Writing expected trait table to:", filename

        f=open(filename,"w")
        f.write("".join(exp_trait_table_lines))
        f.close()

        #Output a transposed, BIOM format expectation table for comparison with predict_traits output

        #NOTE: this is a clumsy way of getting the translated trait table
        # but more elegant, direct methods (directly feeding data to biom's table_factory)
        # weren't working for me readily.   In the future, we should streamline this process
        # Leaving as is for now since this code is mostly for developers so speed/elegence
        # are probably not essential here.

        #Let the hackishness begin

        #Reload the tab-delimited trait table
        header, fields = parse_trait_table(open(filename,"U"))
        fields = [f for f in fields] #converts generator to list

        #Transpose table for .BIOM format so that Observation ids are KOs
        transposed_header, transposed_trait_table_lines =\
          transpose_trait_table_fields(fields,header,\
          id_row_idx=0, input_header_delimiter="\t",output_delimiter="\t")

        #Eliminate newline in header
        trans_trait_table_lines = [transposed_header.strip()]
        trans_trait_table_lines.extend(["\t".join(r) for r in transposed_trait_table_lines])
        trans_trait_table = '\n'.join(trans_trait_table_lines)

        #Write BIOM format expected trait table
        base_name = "--".join(map(str,["exp_biom_traits",opts.method,curr_dist,safe_tip_to_predict]))

        expected_biom_table = parse_table_to_biom(trans_trait_table.split('\n'),\
            table_format = "tab-delimited")

        #print "Expected_trait_table_lines:",exp_trait_table_lines
        filename=os.path.join(opts.output_dir,base_name)
        if opts.verbose:
            print "Writing BIOM-format expected trait table to:", filename

        write_biom_table(expected_biom_table, filename)


        #Write test trait table
        test_trait_table_fields = test_trait_table_fields
        if expected_traits in test_trait_table_fields:
            test_trait_table_fields.remove(expected_traits)
        test_trait_table_lines = [trait_table_header]
        test_trait_table_lines.extend(["\t".join(r)+"\n" for r in test_trait_table_fields])

        #print "Test_trait_table_lines:",test_trait_table_lines
        base_name = "--".join(map(str,["test_trait_table",opts.method,curr_dist,safe_tip_to_predict]))
        filename=os.path.join(opts.output_dir,base_name)

        if opts.verbose:
            print "Writing test trait table to:", filename

        f=open(filename,"w")
        f.write("".join(test_trait_table_lines))
        f.close()

    if opts.verbose:
        print "Done generating test datasets"




    #    print "\n"
    #    print "Original Tree:\n", tree.asciiArt()
    #    print "\n"
    #    print "Test Tree:\n",test_tree.asciiArt()


if __name__ == "__main__":
    main()
