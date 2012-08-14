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
        make_option('-p','--pool_by',\
          default='distance',help='pass comma-separated categories to pool results by those metadata categories. Valid categories are: holdout_method, prediction_method,weighting_method,distance and organism. For example, pass "distance" to output results pooled by holdout distance in addition to holdout method and prediction method  [default: %default]')
]
script_info['version'] = __version__


def unzip(l):
    """Unzip a list into a tuple"""
    return tuple(zip(*l))

def evaluate_test_dataset_dir(obs_dir_fp,exp_dir_fp,file_name_delimiter="--",\
        file_name_field_order=\
        {'file_type':0,"prediction_method":1,"weighting_method":2,"holdout_method":3,\
          "distance":4,"organism":5},strict=False, verbose=True,pool_by=['distance'],\
          roc_success_criteria=['binary','exact']):
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

    pool_by -- if passed, concatenate traits from each trial that is identical in this category.  e.g. pool_by 'distance' will pool traits across individual test genomes with the same holdoout distance.
    
    roc_success_criteria -- a list of methods for measuring 'success' of a prediction.  Separate  ROC curves will be created for each.
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
    #TODO: separate out into a 'get_paired_data_from_dirs' function

    pooled_observations = {}
    pooled_expectations = {}
    for f in sorted(listdir(obs_dir_fp)):
        if verbose:
            print "Examining file: %s" %f
        filename_components = f.split(file_name_delimiter)
        try:
            file_type,holdout_method,weighting_method,\
            prediction_method,distance,organism = \
              filename_components[file_name_field_order['file_type']],\
              filename_components[file_name_field_order['holdout_method']],\
              filename_components[file_name_field_order['weighting_method']],\
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
        base_tag =  '%s\t%s\t' %(holdout_method,prediction_method)
        tags = [base_tag+'all_results']
        if pool_by: 
            combined_tag = base_tag +\
                    "\t".join([str(field)+"_"+str(filename_components[file_name_field_order[field]]) for field in pool_by])
            tags.append(combined_tag)
        
        if verbose:
          print "Pooling by:", pool_by
          print "Combined tags:",tags
        

        #Update the global pooled obs,exp tables with the current results
        #Format results for printing
        for tag in tags:
            #pooled entries will either be empty or be a .BIOM table
            if pooled_observations.get(tag,None) is None:
                pooled_observations[tag] = obs_table
                pooled_expectations[tag] = obs_table
            else:    
                #Entry should already be a table.  So we want to update it by merging in 
                #the current result
                if verbose:
                
                    print "Merging observations with existing biom table for ", tag
                old_obs_table = pooled_observations[tag]
                pooled_observations[tag] =\
                    old_obs_table.merge(obs_table,Sample='union',Observation='union')
                
                if verbose:
                    print "Merging observations with existing biom table for ", combined_tag
                
                old_exp_table = pooled_expectations[tag]
                pooled_expectations[tag] =\
                    old_exp_table.merge(exp_table,Sample='union',Observation='union')

    return run_accuracy_calculations_on_pooled_data(pooled_observations,\
      pooled_expectations,roc_success_criteria=roc_success_criteria,verbose=verbose)


def run_accuracy_calculations_on_pooled_data(pooled_observations,pooled_expectations,roc_success_criteria=['binary','exact'],verbose=False):
    """Run pearson, spearman calculations on pooled observation & expectation datai
    
    pooled_observations -- a dict of observations, with keys set to a descriptive tag and values equal to the observed .biom Table object

    pooled_expectations -- a dict of expectations, with keys set to a descriptive tag and values equal to the observed .biom Table object
    success_criterion -- criterion for success in ROC trial.  Either 'binary' or 'exact'
    verbose -- if set to True, print verbose output
    """
    all_scatter_lines = []
    all_correlation_lines = []
    trials = defaultdict(list)
    for tag in sorted(pooled_observations.keys()):
        if verbose:
            print "calculating scatter,correlations,trials for tag:",tag
        
        metadata= tag.split('\t')
        obs_table = pooled_observations[tag]
        exp_table = pooled_expectations[tag]
        scatter_lines,correlation_lines,new_trial =\
           run_accuracy_calculations_on_biom_table(obs_table,\
           exp_table,metadata,verbose=verbose)
        #trials["_".join(map(str,[holdout_method,prediction_method]))].append(new_trial)
        trials["_".join(map(str,[metadata[0],metadata[1]]))].append(new_trial)
        #field[0] = holdout_method, field[1] = prediction_method
        all_scatter_lines.extend(scatter_lines)
        all_correlation_lines.extend(correlation_lines)
    
    if verbose:
        print "Running ROC analysis..."
    
    
    # Now that we have all trials calculated, we can produce AUC results
    # (ROC stands for receiver-operating characteristics)

    #TODO: calculate the 'binary' and 'exact' ROC curves for each dataset
    all_roc_result_lines={}
    all_roc_auc_lines={}
    for success_criterion in roc_success_criteria:
        if verbose:
            print "Calculating ROC graph points, AUC for criterion:",success_criterion
        roc_result_lines, roc_auc_lines = run_and_format_roc_analysis(trials,success_criterion=success_criterion,verbose=verbose)
        if verbose:
            for l in roc_auc_lines:
                print l
        all_roc_result_lines[success_criterion]=roc_result_lines
        all_roc_auc_lines[success_criterion]=roc_auc_lines
        
    return all_scatter_lines, all_correlation_lines, all_roc_result_lines, all_roc_auc_lines

  
def run_accuracy_calculations_on_biom_table(obs_table,exp_table,metadata,verbose=False): 
        

    scatter_data_points, correlations=\
        evaluate_test_dataset(obs_table,exp_table)
        
    #For AUC, format = [(all_obs_points,all_exp_points)]
    new_trial = unzip(scatter_data_points)
        
    new_scatter_lines = format_scatter_data(scatter_data_points,metadata)
    #scatter_lines.extend(new_scatter_lines)

    new_correlation_lines = format_correlation_data(correlations,metadata)
    #correlation_lines.extend(new_correlation_lines)

    if verbose:
        for l in new_correlation_lines:
            print l
    return new_scatter_lines,new_correlation_lines,new_trial 

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    pool_by = opts.pool_by.split(',') 
    
    file_name_field_order={'file_type':0,"prediction_method":1,\
      "weighting_method":2,"holdout_method":3,"distance":4,"organism":5}

    for k in pool_by:
        
        if k not in file_name_field_order.keys():
            err_text=\
              "Bad value for option '--pool_by'.  Can't pool by '%s'.   Valid categories are: %s" %(k,\
              ",".join(file_name_field_order.keys()))
            raise ValueError(err_text)
    if opts.verbose:
        print "Pooling results by:",pool_by
    
    
    roc_success_criteria = ['binary','exact','int_exact']

    scatter_lines,correlation_lines,roc_result_lines,roc_auc_lines =\
      evaluate_test_dataset_dir(opts.trait_table_dir,\
      opts.exp_trait_table_dir,file_name_delimiter="--",\
      file_name_field_order=file_name_field_order,pool_by=pool_by,\
      roc_success_criteria=roc_success_criteria,verbose=opts.verbose)
   

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
    for c in roc_result_lines.keys(): 
        output_fp = join(opts.output_dir,'evaluation_roc_data_%s.tab' %c)
        if opts.verbose:
            print "Outputting ROC data for success criterion %s to: %s" %(c,output_fp)
        file_lines = roc_result_lines[c]
    
        f = open(output_fp,"w+")
        f.writelines(file_lines)
        f.close()

    #Output summary ROC AUC data
    if opts.verbose:
        print "Writing ROC AUC data..."
    
    for c in roc_auc_lines.keys(): 
        output_fp = join(opts.output_dir,'evaluation_roc_auc_data_%s.tab' %c)
        file_lines = roc_auc_lines[c]
    
        if opts.verbose:
            print "Outputting ROC AUC data for success criterion %s to: %s" %(c,output_fp)
        f = open(output_fp,"w+")
        f.writelines(file_lines)
        f.close()



if __name__ == "__main__":
    main()
