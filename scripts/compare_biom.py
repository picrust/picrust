#!/usr/bin/env python
from __future__ import division

__author__ = "Morgan Langille"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "1.1.3"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"

from collections import defaultdict
from os import listdir
from numpy import array,ravel,asarray
from os.path import join,basename
from cogent.util.option_parsing import parse_command_line_parameters,\
  make_option
from picrust.evaluate_test_datasets import calculate_accuracy_stats_from_observations
from biom import load_table, Table
from picrust.util import make_output_dir_for_file
from random import shuffle

script_info = {}
script_info['brief_description'] = "Compare the accuracy of biom files (expected and observed) either by observations (default) or by samples."
script_info['script_description'] =\
    """ """
script_info['script_usage'] = [\
    ("Example 1","Compare an observed table to an expected table using relative abundance","%prog -e expected_ra.biom -o compare_results_ra.tab observed_ra.biom"),
    ("Example 2","Compare an observed table to an expected table using real counts","%prog --not_relative_abundance -e expected.biom -o compare_results.tab observed.biom")]
script_info['output_description']= "Outputs will be tab delimited file with various accuracy metrics."
script_info['required_options'] = [
 make_option('-e','--exp_trait_table_fp',type="existing_filepath",help='the expected trait table (biom format)'),\
 make_option('-o','--output_fp',type="new_filepath",help='the output file'),
]
script_info['optional_options'] = [
  make_option('-c','--compare_observations',action="store_true",default=False,help='Calculate accuracy values by comparing between observations (instead of between samples) [default: %default]'),\
  make_option('-n','--normalize',action="store_true",default=False,help='Convert both expected and observed tables to relative abundances (instead of observations) [default: %default]'),
  make_option('-l','--limit_to_expected_observations',action="store_true",default=False,help='Ignore observations that are not in the expected table[default: %default]'),
  make_option('--limit_to_observed_observations',action="store_true",default=False,help='Ignore observations that are not in the observed table[default: %default]'),
  make_option('-s','--shuffle_samples',action="store_true",default=False,help='Shuffle samples ids randomly before measuring accuracy[default: %default]'),
  make_option('--not_relative_abundance_scores',action="store_true",default=False,help='Round numbers (instead of taking ceil() which is used for RA) before calculating TP,FP,FN,TN [default: %default]')

        ]
script_info['disallow_positional_arguments'] = False
script_info['version'] = __version__


def match_biom_tables(observed_table,expected_table_keep,verbose=False,limit_to_expected_observations=False,limit_to_observed_observations=False,normalize=False,shuffle_samples=False):

    expected_table = expected_table_keep.copy()

    overlapping_obs_ids = list(set(observed_table.ids(axis='observation')) &
                               set(expected_table.ids(axis='observation')))

    if len(overlapping_obs_ids) < 1:
        print "obs ids:", observed_table.ids(axis='observation')[0:10]
        print "exp ids:", expected_table.ids(axis='observation')[0:10]

        raise ValueError,\
         "No observation ids are in common  between the observed and expected tables, so no evaluations can be performed."


    if limit_to_expected_observations:
        def f(data_vector, id_, metadata):
            return (id_ in overlapping_obs_ids)
        observed_table=observed_table.filter(f, axis='observation', inplace=False)

    if limit_to_observed_observations:
        def f(data_vector, id_, metadata):
            return (id_ in overlapping_obs_ids)
        expected_table=expected_table.filter(f, axis='observation', inplace=False)


    ###Make tables have same set (e.g.number) of ObservationIds and in the same order###
    #1)identify ObservationIds unique to each table
    unique_obs_in_expected=list(set(expected_table.ids(axis='observation')) - set(observed_table.ids(axis='observation')))
    unique_obs_in_observed=list(set(observed_table.ids(axis='observation')) - set(expected_table.ids(axis='observation')))

    #2)Add each missing observation with all 0's

    if unique_obs_in_observed:
        empty_obs_data=[[0]*len(expected_table.ids())]*len(unique_obs_in_observed)
        empty_obs_data=asarray(empty_obs_data)
        empty_obs_table=Table(empty_obs_data, unique_obs_in_observed,
                              expected_table.ids())
        expected_table=expected_table.merge(empty_obs_table)

    if unique_obs_in_expected:
        empty_obs_data=[[0]*len(observed_table.ids())]*len(unique_obs_in_expected)
        empty_obs_data=asarray(empty_obs_data)
        empty_obs_table=Table(empty_obs_data, unique_obs_in_expected,
                              observed_table.ids())
        observed_table=observed_table.merge(empty_obs_table)


    #3)sort the ObservationIds so they are in the same order between the tables

    if verbose:
        print "Sorting observations in expected table to match observed table..."
    expected_table=expected_table.sort_order(observed_table.ids(axis='observation'),
                                             axis='observation')


    overlapping_sample_ids = list(set(observed_table.ids()) &
                            set(expected_table.ids()))

    if verbose:
        num_uniq_obs_sample_ids=len(observed_table.ids())-len(overlapping_sample_ids)
        num_uniq_exp_sample_ids=len(expected_table.ids())-len(overlapping_sample_ids)
        if num_uniq_obs_sample_ids:
            print "Num observed samples not in expected: {0}".format(num_uniq_obs_sample_ids)
        if num_uniq_exp_sample_ids:
            print "Num expected samples not in observed: {0}".format(num_uniq_exp_sample_ids)
        print "Num samples with same id: {0}".format(len(overlapping_sample_ids))


    if normalize:
        if verbose:
            print "Normalizing tables..."
        observed_table=observed_table.norm(axis='sample', inplace=False)
        expected_table=expected_table.norm(axis='sample', inplace=False)

    if verbose:
        print "Extracting data from biom objects..."
    # create lists to contain filtered data - we're going to need the data in
    # numpy arrays, so it makes sense to compute this way rather than filtering
    # the tables
    obs_data = {}
    exp_data = {}

    # build lists of filtered data
    for sample_id in overlapping_sample_ids:
        exp_data[sample_id]=expected_table.data(sample_id)


    if shuffle_samples:
        if verbose:
            print "Randomly shufflying sample ids..."
        sample_ids_to_shuffle=overlapping_sample_ids[:]
        shuffle(sample_ids_to_shuffle)

        for index in range(len(overlapping_sample_ids)):
            obs_data[overlapping_sample_ids[index]]=observed_table.data(sample_ids_to_shuffle[index])
    else:
        for sample_id in overlapping_sample_ids:
            obs_data[sample_id]=observed_table.data(sample_id)

    return obs_data,exp_data



def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    verbose=opts.verbose

    min_args = 1
    if len(args) < min_args:
       option_parser.error('One or more predicted biom files must be provided.')
    observed_files=args


    make_output_dir_for_file(opts.output_fp)
    out_fh=open(opts.output_fp,'w')

    if verbose:
        print "Loading expected trait table file:",opts.exp_trait_table_fp

    exp_table = load_table(opts.exp_trait_table_fp)

    header_printed=False
    header_keys=[]
    delimiter="\t"


    for observed_file in observed_files:
        observed_file_name = basename(observed_file)

        if verbose:
            print "Loading predicted trait table file:",observed_file_name

        obs_table = load_table(observed_file)

        if opts.compare_observations:
            if verbose:
                print "Transposing tables to allow evaluation of observations (instead of samples)..."
            obs_table = obs_table.transpose()
            exp_table = exp_table.transpose()

        if verbose:
           print "Matching predicted and expected tables..."

        obs,exp=match_biom_tables(obs_table,exp_table,verbose=verbose,limit_to_expected_observations=opts.limit_to_expected_observations,limit_to_observed_observations=opts.limit_to_observed_observations,normalize=opts.normalize,shuffle_samples=opts.shuffle_samples)

        if verbose:
            print "Calculating accuracy stats for all observations..."

        #import pdb; pdb.set_trace()
        for i in obs:
            if verbose:
                print "Calculating stats for: ",i
            if opts.not_relative_abundance_scores:
                results=calculate_accuracy_stats_from_observations(obs[i],exp[i],success_criterion='binary')
            else:
                results=calculate_accuracy_stats_from_observations(obs[i],exp[i],success_criterion='ra_exact')

            #If first pass then print out header
            if not header_printed:
                header_printed=True
                header_keys=sorted(results.keys())
                out_fh.write(delimiter.join(['file','label']+header_keys)+"\n")

            #print results using same order as header
            values=[observed_file_name,i]+['{0:.3g}'.format(results[x]) for x in header_keys]
            out_str=delimiter.join(map(str,values))+"\n"
            out_fh.write(out_str)


if __name__ == "__main__":
    main()
