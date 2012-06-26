#!/usr/bin/env python
# File created on 10 April 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2012, The PICRUST project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"
 
from collections import defaultdict
from os import listdir
from os.path import join
from cogent.util.option_parsing import parse_command_line_parameters,\
  make_option
from picrust.evaluate_test_datasets import evaluate_test_dataset,\
  format_scatter_data, format_correlation_data, run_and_format_roc_analysis

from biom.parse import parse_biom_table, convert_biom_to_table


from picrust.parse import parse_trait_table, parse_marker_gene_copy_numbers

script_info = {}
script_info['brief_description'] = "Evaluate the accuracy of character predictions, given directories of expected vs. observed test results"
script_info['script_description'] =\
    """The script finds all paired expected and observed values in a set of directories and generates the following output: 1) data for a scatterplot of observed vs. expected values for each character (typically gene family count) within each organism (so one file per organism). 2) A summary of accuracy across all organisms.   
    character """
script_info['script_usage'] = [("","Evaluate the accuracy of all predictions in a folder, and output summary statistics.","%prog -i obs_otu_table.biom -e exp_otu_table.txt -o./evaluation_results/")]
script_info['output_description']= "Outputs will be obs,exp data points for the comparison"
script_info['required_options'] = [
 make_option('-i','--trait_table_dir',type="existing_dirpath",help='the input trait table directory (files in biom format)'),\
 make_option('-e','--exp_trait_table_dir',type="existing_dirpath",help='the input expected trait table directory (files in biom format)'),\
 make_option('-o','--output_dir',type="new_dirpath",help='the output directory'),
]
script_info['optional_options'] = [
        #make_option('--test_trait_table_format',choices=['biom','tab-delimited'],\
        #  default='tab-delimited',help='the format of the test trait tables. Choices are: %choices [default: %default]'),\
        #make_option('--expected_trait_table_format',choices=['biom','tab-delimited'],\
        #  default='tab-delimited',help='the format of the expected trait tables. Choices are: %choices [default: %default]')

]
script_info['version'] = __version__



def evaluate_test_dataset_dir(obs_dir_fp,exp_dir_fp,file_name_delimiter="--",\
        file_name_field_order=\
          {'file_type':0,"prediction_method":1, "holdout_method":2,\
          "distance":3,"organism":4},strict=False, verbose=True):
    """Return control evaluation results from the given directories
    
    obs_dir_fp -- directory containing PICRUST-predicted genomes.   These MUST start with
    'predict_traits', and must contain the values specified in file_name_field_order,\
    separated by the delimiter given in file_name_delimiter.  For example:

    predict_traits--exclude_tips_by_distance--0.87--'NC_000913|646311926'

    exp_dir_fp -- as obs_dir_fp above, but expectation file names (usually sequenced genomes
    with known gene content) must start with exp_biom_traits

    file_name_delimiter -- the delimiter that separates metadata stored in the filename
    
    NOTE: technically this isn't the best way of doing things.  We may want at some point
    to revisit this setup and store metadata about each comparison in a separate file.  But
    storing in the filename is convenient for our initial analysis.

    file_name_field_order -- the order of the required metadata fields in the filename.
    Required fields are file_type,method,distance,and organism


    Description:

    The method assumes that for each file type in the observed directory, a paired file
    is also found in the exp_dir with similar method, distance, and organism, but a varied 
    file type (test_tree, test_trait_table)

    
    Process:
    1. Search test directory for all gene predictions in the correct format
    2. For each, find the corresponding expected trait table in the expectation file
    3. Organize all of the correlations and scatter points
    4. Return the following outputs:
        --- Table of obs,exp values with organism, method, distance metadata
        --- Summary of correlations (table by organism and distance)
        --- Summary of AUC (table by organism,method,distance)


    
    
    """
    
    
    trials = defaultdict(list)
    correlation_lines = []
    scatter_lines = []

    #We'll want a quick unzip fn for converting points to trials
    unzip = lambda l:tuple(zip(*l))
    
    for f in sorted(listdir(obs_dir_fp)):
        if verbose:
            print "Examining file: %s" %f
        filename_components = f.split(file_name_delimiter)
        try:
            file_type,holdout_method,prediction_method,distance,organism = \
              filename_components[file_name_field_order['file_type']],\
              filename_components[file_name_field_order['holdout_method']],\
              filename_components[file_name_field_order['prediction_method']],\
              filename_components[file_name_field_order['distance']],\
              filename_components[file_name_field_order['organism']]
        except IndexError, e:
            print "Could not parse filename %s using delimiter: %s.  Skipping..." %(f,file_name_delimiter)
            continue
        if verbose:
            print "HOLDOUT METHOD:", holdout_method

        #Get predicted traits
        if file_type == 'predict_traits':
            if verbose:
                print "Found a prediction file"
                print "Loading .biom format observation table:",f
            
            #try:
            obs_table =\
              parse_biom_table(open(join(obs_dir_fp,f),'U'))
            #except:
            #    raise RuntimeError(\
            #      "Could not parse predicted trait file: %s.   Is it a .biom formatted file?" %(f))
        else:
            continue
        
        
        # Get paired observation file
        exp_filename = file_name_delimiter.join(['exp_biom_traits',holdout_method,distance,organism])
        exp_filepath = join(exp_dir_fp,exp_filename)
        if verbose:
            print "Looking for the expected trait file matching %s here: %s" %(f,exp_filepath)

        try:
            exp_table =\
              parse_biom_table(open(exp_filepath,"U"))
        except IOError, e:
            if strict:
                raise IOError(e)
            else:
                if verbose:
                    print "Missing expectation file....skipping!"
                continue

        scatter_data_points, correlations=\
          evaluate_test_dataset(obs_table,exp_table)
        
        #For AUC, format = [(all_obs_points,all_exp_points)]
        new_trial = unzip(scatter_data_points)
        trials["_".join(map(str,[holdout_method,prediction_method]))].append(new_trial)
        #Format results for printing
        metadata = [organism,holdout_method,prediction_method,distance]
        
        new_scatter_lines = format_scatter_data(scatter_data_points,metadata)
        scatter_lines.extend(new_scatter_lines)

        new_correlation_lines = format_correlation_data(correlations,metadata)
        correlation_lines.extend(new_correlation_lines)

        if verbose:
            for l in new_correlation_lines:
                print l
        
    # Now that we have all trials calculated, we can produce AUC results
    # (ROC stands for receiver-operating characteristics)

    roc_result_lines, roc_auc_lines = run_and_format_roc_analysis(trials)
    if verbose:
        for l in roc_auc_lines:
            print l
    return scatter_lines, correlation_lines, roc_result_lines, roc_auc_lines

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    scatter_lines,correlation_lines,roc_result_lines,roc_auc_lines =\
      evaluate_test_dataset_dir(opts.trait_table_dir,\
      opts.exp_trait_table_dir,file_name_delimiter="--",\
      file_name_field_order={'file_type':0,"prediction_method":1,\
      "holdout_method":2,"distance":3,"organism":4})
   

    #Output scatter data
    if opts.verbose:
        print "Writing scatter plot data..."
    
    output_fp = join(opts.output_dir,'evaluation_scatter_data.tab')
    file_lines = scatter_lines
    
    f = open(output_fp,"w+")
    f.writelines(file_lines)
    f.close()

    #Output correlation data
    if opts.verbose:
        print "Writing correlation data..."
    
    output_fp = join(opts.output_dir,'evaluation_correlation_data.tab')
    file_lines = correlation_lines
    
    f = open(output_fp,"w+")
    f.writelines(file_lines)
    f.close()

    #Output raw ROC plot data
    if opts.verbose:
        print "Writing ROC data..."
    
    output_fp = join(opts.output_dir,'evaluation_roc_data.tab')
    file_lines = roc_result_lines
    
    f = open(output_fp,"w+")
    f.writelines(file_lines)
    f.close()

    #Output summary ROC AUC data
    if opts.verbose:
        print "Writing ROC AUC data..."
    
    output_fp = join(opts.output_dir,'evaluation_roc_auc_data.tab')
    file_lines = roc_auc_lines
    
    f = open(output_fp,"w+")
    f.writelines(file_lines)
    f.close()



if __name__ == "__main__":
    main()
