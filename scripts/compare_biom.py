#!/usr/bin/env python
from __future__ import division

__author__ = "Morgan Langille"
__copyright__ = "Copyright 2012, The PICRUST project"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"
 
from collections import defaultdict
from os import listdir
from numpy import array,ravel
from os.path import join,basename
from cogent.util.option_parsing import parse_command_line_parameters,\
  make_option
from picrust.evaluate_test_datasets import calculate_accuracy_stats_from_observations
from biom.parse import parse_biom_table, convert_biom_to_table
from picrust.util import make_output_dir

from biom.table import table_factory,DenseOTUTable,SparseOTUTable
from picrust.parse import parse_trait_table, parse_marker_gene_copy_numbers

script_info = {}
script_info['brief_description'] = "Compare the accuracy of biom files (expected and observed) either by observations (default) or by samples."
script_info['script_description'] =\
    """ """
script_info['script_usage'] = [("","","%prog -e exp_otu_table.biom -o results.tab obs_otu_table.biom [obs_otu_table2.biom]")]
script_info['output_description']= "Outputs will be tab delimited file with various accuracy metrics."
script_info['required_options'] = [
 make_option('-e','--exp_trait_table_fp',type="existing_filepath",help='the expected trait table (biom format)'),\
 make_option('-o','--output_fp',type="new_filepath",help='the output file'),
]
script_info['optional_options'] = [
  make_option('-c','--compare_observations',action="store_true",default=False,help='Calculate accuracy values by comparing between observations (instead of between samples) [default: %default]'),\
  make_option('-n','--normalize',action="store_true",default=False,help='Convert both expected and observed tables to relative abundances (instead of observations) [default: %default]'),
  make_option('-l','--limit_to_expected_observations',action="store_true",default=False,help='Ignore observations that are not in the expected table[default: %default]')
  
        ]
script_info['disallow_positional_arguments'] = False
script_info['version'] = __version__




def match_biom_tables_by_observation(observed_table,expected_table_keep,verbose=False):

    expected_table = expected_table_keep.copy()
    overlapping_obs_ids = list(set(observed_table.ObservationIds) &
                            set(expected_table.ObservationIds))
   
    if len(overlapping_obs_ids) < 1:
        print "obs ids:",observed_table.ObservationIds[0:10]
        print "exp ids:",expected_table.ObservationIds[0:10]
        
        raise ValueError,\
         "No observation ids are in common  between the observed and expected tables, so no evaluations can be performed."

    overlapping_sample_ids = list(set(observed_table.SampleIds) &
                            set(expected_table.SampleIds))


    if len(overlapping_sample_ids) != len(observed_table.SampleIds) or len(overlapping_sample_ids) != len(expected_table.SampleIds):
        print "obs ids:",observed_table.SampleIds
        print "exp ids:",expected_table.SampleIds
        
        raise ValueError,\
         "Tables do not have the same samples. Can't currently compare these types of tables."

    #ensure samples are in same order
    if verbose:
        print "Sorting samples in predicted table..."
    observed_table=observed_table.sortSampleOrder(overlapping_sample_ids)

    if verbose:
        print "Sorting samples in expected table..."
    expected_table=expected_table.sortSampleOrder(overlapping_sample_ids)

    # create lists to contain filtered data - we're going to need the data in 
    # numpy arrays, so it makes sense to compute this way rather than filtering
    # the tables
    obs_data = {}
    exp_data = {}
    
    #observed_table=observed_table.normObservationBySample()
    
    # build lists of filtered data
    for obs_id in overlapping_obs_ids:
        obs_data[obs_id]=observed_table.observationData(obs_id)
        exp_data[obs_id]=expected_table.observationData(obs_id)

    #print obs_data
    #print exp_data
#    flat_obs_data = ravel(array(obs_data))
#    flat_exp_data = ravel(array(exp_data))
  #  return flat_obs_data,flat_exp_data
    return obs_data,exp_data

def match_biom_tables_by_sample(observed_table,expected_table_keep,verbose=False,limit_to_expected_observations=False,normalize=False):
   
    expected_table = expected_table_keep.copy()
    overlapping_obs_ids = list(set(observed_table.ObservationIds) &
                            set(expected_table.ObservationIds))
   
    if len(overlapping_obs_ids) < 1:
        print "obs ids:",observed_table.ObservationIds[0:10]
        print "exp ids:",expected_table.ObservationIds[0:10]
        
        raise ValueError,\
         "No observation ids are in common  between the observed and expected tables, so no evaluations can be performed."

    #import pdb; pdb.set_trace()

    if limit_to_expected_observations:
        def f(data_vector, id_, metadata):
            return (id_ in overlapping_obs_ids)
        observed_table=observed_table.filterObservations(f)
        

    ###Make tables have same set (e.g.number) of ObservationIds and in the same order###
    #1)identify ObservationIds unique to each table
    unique_obs_in_expected=list(set(expected_table.ObservationIds) - set(observed_table.ObservationIds))
    unique_obs_in_observed=list(set(observed_table.ObservationIds) - set(expected_table.ObservationIds))
    
    #2)Add each missing observation with all 0's

    if unique_obs_in_observed:
        empty_obs_data=[[0]*len(expected_table.SampleIds)]*len(unique_obs_in_observed)
        empty_obs_table=table_factory(empty_obs_data,expected_table.SampleIds,unique_obs_in_observed, constructor=DenseOTUTable)
        expected_table=expected_table.merge(empty_obs_table)

    if unique_obs_in_expected:
        empty_obs_data=[[0]*len(observed_table.SampleIds)]*len(unique_obs_in_expected)
        empty_obs_table=table_factory(empty_obs_data,observed_table.SampleIds,unique_obs_in_expected, constructor=DenseOTUTable)
        observed_table=observed_table.merge(empty_obs_table)
       

    #3)sort the ObservationIds so they are in the same order between the tables
   
    if verbose:
        print "Sorting observations in expected table to match observed table..."
    expected_table=expected_table.sortObservationOrder(observed_table.ObservationIds)


    overlapping_sample_ids = list(set(observed_table.SampleIds) &
                            set(expected_table.SampleIds))


    if normalize:
        if verbose:
            print "Normalizing tables..."
        observed_table=observed_table.normObservationBySample()        
        expected_table=expected_table.normObservationBySample()        


    # create lists to contain filtered data - we're going to need the data in 
    # numpy arrays, so it makes sense to compute this way rather than filtering
    # the tables
    obs_data = {}
    exp_data = {}
    
    # build lists of filtered data
    for sample_id in overlapping_sample_ids:
        obs_data[sample_id]=observed_table.sampleData(sample_id)
        exp_data[sample_id]=expected_table.sampleData(sample_id)

    return obs_data,exp_data
        

   
def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    verbose=opts.verbose

    min_args = 1
    if len(args) < min_args:
       option_parser.error('One or more predicted biom files must be provided.')
    observed_files=args
   

    if verbose:
        print "Loading expected trait table file:",opts.exp_trait_table_fp

    exp_table =parse_biom_table(open(opts.exp_trait_table_fp,'U'))


    out_fh=open(opts.output_fp,'w')

    header_printed=False
    header_keys=[]
    delimiter="\t"


    for observed_file in observed_files:
        observed_file_name=basename(observed_file)

        if verbose:
            print "Loading predicted trait table file:",observed_file_name

        obs_table =parse_biom_table(open(observed_file,'U'))

        if verbose:
           print "Matching predicted and expected tables..."    

        if opts.compare_observations:
            obs,exp=match_biom_tables_by_observation(obs_table,exp_table,verbose=verbose)
        else:
            obs,exp=match_biom_tables_by_sample(obs_table,exp_table,verbose=verbose,limit_to_expected_observations=opts.limit_to_expected_observations,normalize=opts.normalize)
           
        if verbose:
            print "Calculating accuracy stats for all observations..."


        for i in obs:
            if verbose:
                print "Calculating stats for: ",i
            results=calculate_accuracy_stats_from_observations(obs[i],exp[i],success_criterion='ra_exact')

            #If first pass then print out header
            if not header_printed:
                header_printed=True
                header_keys=sorted(results.keys())
                out_fh.write(delimiter.join(['file','label']+header_keys)+"\n")

            #print results using same order as header
            values=[observed_file_name,i]+[results[x] for x in header_keys]
            out_str=delimiter.join(map(str,values))+"\n"
            out_fh.write(out_str)


if __name__ == "__main__":
    main()
