#!/usr/bin/env python
# File created on 10 April 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
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
from picrust.evaluate_test_datasets import unzip,evaluate_test_dataset,\
 update_pooled_data, run_accuracy_calculations_on_biom_table,run_accuracy_calculations_on_pooled_data,\
 format_scatter_data, format_correlation_data, run_and_format_roc_analysis

from biom.parse import parse_biom_table, convert_biom_to_table


from picrust.parse import parse_trait_table, parse_marker_gene_copy_numbers

script_info = {}
script_info['brief_description'] = "Pool character predictions within a directory, given directories of expected vs. observed test results"
script_info['script_description'] =\
    """The script finds all paired expected and observed values in a set of directories and generates pooled .biom files in a specified output directory"""
script_info['script_usage'] = [("","Pool .biom files according to holdout_distance.","%prog -i obs_otu_table_dir -e exp_otu_table_dir -p distance -o./evaluation_results/pooled_by_distance/")]
script_info['output_description']= "Outputs will be obs,exp data points for the comparison"
script_info['required_options'] = [
 make_option('-i','--trait_table_dir',type="existing_dirpath",help='the input trait table directory (files in biom format)'),\
 make_option('-e','--exp_trait_table_dir',type="existing_dirpath",help='the input expected trait table directory (files in biom format)'),\
 make_option('-o','--output_dir',type="new_dirpath",help='the output directory'),
]
script_info['optional_options'] = [
        make_option('-f','--field_order',\
                default='file_type,prediction_method,weighting_method,holdout_method,distance,organism',help='pass comma-separated categories, in the order they appear in file names.   Categories are "file_type","prediction_method","weighting_method","holdout_method" (randomization vs. holdout),"distance",and "organism".  Example:  "-f file_type,test_method,asr_method specifies that files will be in the form: predict_traits--distance_exclusion--wagner.  Any unspecified values are set to "not_specified".  [default: %default]'),\
        make_option('-p','--pool_by',\
          default=False,help='pass comma-separated categories to pool results by those metadata categories. Valid categories are: holdout_method, prediction_method,weighting_method,distance and organism. For example, pass "distance" to output results pooled by holdout distance in addition to holdout method and prediction method  [default: %default]')
]
script_info['version'] = __version__


def iter_prediction_expectation_pairs(obs_dir_fp,exp_dir_fp,file_name_field_order,file_name_delimiter,verbose=False):
    """Iterate pairs of observed, expected biom file names"""
    input_files=sorted(listdir(obs_dir_fp))
    for file_number,f in enumerate(input_files):
        if verbose:
            print "\nExamining file {0} of {1}: {2}".format(file_number+1,len(input_files),f)
        if 'accuracy_metrics' in f:
            print "%s is an Accuracy file...skipping" %str(f)
            continue
        #filename_components_list = f.split(file_name_delimiter)
        #Get predicted traits
        filename_metadata = get_metadata_from_filename(f,file_name_field_order,\
          file_name_delimiter,verbose=verbose)

        if filename_metadata.get('file_type',None) == 'predict_traits':
            if verbose:
                #print "Found a prediction file"
                print "\tLoading .biom format observation table:",f
            
            try:
              obs_table =\
                parse_biom_table(open(join(obs_dir_fp,f),'U'))
            except ValueError:
                print 'Failed, skipping...'
                continue
            #    raise RuntimeError(\
            #      "Could not parse predicted trait file: %s.   Is it a .biom formatted file?" %(f))
        else:
            continue
        
        # Get paired observation file
        exp_filename = file_name_delimiter.join(['exp_biom_traits',filename_metadata['holdout_method'],filename_metadata['distance'],filename_metadata['organism']])
        exp_filepath = join(exp_dir_fp,exp_filename)
        if verbose:
            print "\tLooking for the expected trait file matching %s here: %s" %(f,exp_filepath)

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
        yield obs_table,exp_table,f

def get_metadata_from_filename(f,file_name_field_order,file_name_delimiter,\
  default_text='not_specified',verbose=False):
    """Extract metadata values from a filename"""    
    filename_components = {}
    for i,field in enumerate(f.split(file_name_delimiter)):
        filename_components[i]=field
    #if verbose:
    #    print "Filename components:",filename_components
    filename_metadata = {}
    try:
        for field in file_name_field_order.keys():
           filename_metadata[field] =\
             filename_components.get(file_name_field_order.get(field,default_text),default_text)

        #if verbose:
        #    print "filename_metadata:",filename_metadata
    except IndexError, e:
        print "Could not parse filename %s using delimiter: %s.  Skipping..." %(f,file_name_delimiter)
        return None

    return filename_metadata

def pool_test_dataset_dir(obs_dir_fp,exp_dir_fp,file_name_delimiter="--",\
        file_name_field_order=\
        {'file_type':0,"prediction_method":1,"weighting_method":2,"holdout_method":3,\
          "distance":4,"organism":5},strict=False, verbose=True,pool_by=['distance']):
    """Retrun pooled control &  evaluation results from the given directories
    
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

    pool_by -- if passed, concatenate traits from each trial that is identical in this category.  e.g. pool_by 'distance' will pool traits across individual test genomes with the same holdout distance.
    
    The method assumes that for each file type in the observed directory, a paired file
    is also found in the exp_dir with similar method, distance, and organism, but a varied 
    file type (test_tree, test_trait_table)

    
    Process:
    1. Search test directory for all gene predictions in the correct format
    2. For each, find the corresponding expected trait table in the expectation file
    3. Pool by specified pool_by values
    4. Return dicts of pooled observation,expectation values
    """
    trials = defaultdict(list)
    #We'll want a quick unzip fn for converting points to trials
    #TODO: separate out into a 'get_paired_data_from_dirs' function

    pooled_observations = {}
    pooled_expectations = {}
    pairs = iter_prediction_expectation_pairs(obs_dir_fp,exp_dir_fp,file_name_field_order,file_name_delimiter,verbose=verbose)
    file_number = 0
    for obs_table,exp_table,filename in pairs:
        #print "analyzing filename:",filename 
        filename_metadata= get_metadata_from_filename(filename,file_name_field_order,\
          file_name_delimiter,verbose=verbose)
        
        #base_tag =  '%s\t%s\t' %(filename_metadata['holdout_method'],filename_metadata['prediction_method'])
        #tags = [base_tag+'all_results']
        if 'file_type' in pool_by:
            pool_by.remove('file_type') #we do this manually at the end
        combined_tag = ['all']*len(file_name_field_order.keys())
        for field in file_name_field_order.keys():
          #print combined_tag
          #print file_name_field_order
          idx = file_name_field_order[field]
          #print idx
          if field in pool_by:
              combined_tag[idx] = filename_metadata[field]
        tags=[file_name_delimiter.join(combined_tag)]
       
        if verbose:
          print "Pooling by:", pool_by
          print "Combined tags:",tags
        
        pooled_observations,pooled_expectations =\
        update_pooled_data(obs_table,exp_table,tags,pooled_observations,\
        pooled_expectations,str(file_number),verbose=verbose)
        file_number += 1

    return pooled_observations,pooled_expectations 


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    pool_by = opts.pool_by.split(',') 
    
    #Construct a dict from user specified field order
    file_name_field_order = {}
    for i,field in enumerate(opts.field_order.split(',')):
        file_name_field_order[field]=i
        if opts.verbose:
            print "Assuming file names are in this order:",file_name_field_order
    for k in pool_by:
        #Check that we're only pooling by values that exist 
        if k not in file_name_field_order.keys():
            err_text=\
            "Bad value for option '--pool_by'.  Can't pool by '%s'.   Valid categories are: %s" %(k,\
            ",".join(file_name_field_order.keys()))
            raise ValueError(err_text)

        if opts.verbose:
            print "Pooling results by:",pool_by
    
    file_name_delimiter='--'
    pooled_observations,pooled_expectations = pool_test_dataset_dir(opts.trait_table_dir,\
      opts.exp_trait_table_dir,file_name_delimiter=file_name_delimiter,\
      file_name_field_order=file_name_field_order,pool_by=pool_by,\
      verbose=opts.verbose)
    
    #prediction_prefix = 'predict_traits'
    #expectation_prefix = 'exp_biom_traits'

    for tag in pooled_observations.keys():
        obs_table = pooled_observations[tag]
        exp_table = pooled_expectations[tag]

        #obs_table_filename = file_name_delimiter.join([prediction_prefix]+[t for t in tag.split()])
        #exp_table_filename = file_name_delimiter.join([expectation_prefix]+[t for t in tag.split()])
        
        obs_table_filename = file_name_delimiter.join(['predict_traits']+[t for t in tag.split()])
        exp_table_filename = file_name_delimiter.join(['exp_biom_table']+[t for t in tag.split()])
        
        obs_outpath = join(opts.output_dir,obs_table_filename)
        exp_outpath = join(opts.output_dir,exp_table_filename)

        print obs_outpath
        print exp_outpath
        
        f=open(obs_outpath,'w')
        f.write(obs_table.delimitedSelf())
        f.close()

        f=open(exp_outpath,'w')
        f.write(exp_table.delimitedSelf())
        f.close()

   
    

if __name__ == "__main__":
    main()
