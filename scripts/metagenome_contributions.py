#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011, The PICRUST project"
__credits__ = ["Greg Caporaso","Jesse Zaneveld","Morgan Langille"]
__license__ = "GPL"
__version__ = "0.9.1-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"
 


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from biom.parse import parse_biom_table
from picrust.predict_metagenomes import predict_metagenomes, calc_nsti
from picrust.metagenome_contributions import partition_metagenome_contributions
from picrust.util import make_output_dir_for_file, get_picrust_project_dir
from os import path
from os.path import join
import gzip

script_info = {}
script_info['brief_description'] = "This script partitions metagenome functional contributions according to function, OTU, and sample, for a given OTU table."
script_info['script_description'] = ""
script_info['script_usage'] = [
("","Partition the predicted contribution to the  metagenomes from each organism in the given OTU table, limited to only K00001, K00002, and K00004.","%prog -i normalized_otus.biom -l K00001,K00002,K00004 -o ko_metagenome_contributions.tab"),
("","Partition the predicted contribution to the  metagenomes from each organism in the given OTU table, limited to only COG0001 and COG0002.","%prog -i normalized_otus.biom -l COG0001,COG0002 -t COG -o cog_metagenome_contributions.tab")
]
script_info['output_description']= "Output is a tab-delimited column indicating OTU contribution to each function."
script_info['required_options'] = [
 make_option('-i','--input_otu_table',type='existing_filepath',help='the input otu table in biom format'),
 make_option('-o','--output_fp',type="new_filepath",help='the output file for the metagenome contributions')
]
type_of_prediction_choices=['KO','COG']

script_info['optional_options'] = [
        make_option('-t','--type_of_prediction',default='KO',type="choice",\
                    choices=type_of_prediction_choices,\
                    help='Type of functional predictions. Valid choices are: '+\
                    ', '.join(type_of_prediction_choices)+\
                    ' [default: %default]'),
        make_option('-c','--input_count_table',default=None,type="existing_filepath",help='Precalculated function predictions on per otu basis in biom format (can be gzipped). Note: using this option overrides --type_of_prediction. [default: %default]'),
 
        make_option('-l','--limit_to_function',default=None,help='If provided, only output predictions for the specified function ids.  Multiple function ids can be passed using comma delimiters.')
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
  
    if opts.limit_to_function:
        limit_to_functions = opts.limit_to_function.split(',')
        if opts.verbose:
            print "Limiting output to only functions:",limit_to_functions
    else:
        limit_to_functions = []

    if opts.verbose:
        print "Loading otu table: ",opts.input_otu_table

    otu_table = parse_biom_table(open(opts.input_otu_table,'U'))

    if(opts.input_count_table is None):
        if(opts.type_of_prediction == 'KO'):
            input_count_table=join(get_picrust_project_dir(),'picrust','data','ko_precalculated.biom.gz')
        elif(opts.type_of_prediction == 'COG'):
            input_count_table=join(get_picrust_project_dir(),'picrust','data','cog_precalculated.biom.gz')
    else:
        input_count_table=opts.input_count_table

    if opts.verbose:
        print "Loading trait table: ", input_count_table

    ext=path.splitext(input_count_table)[1]

    if opts.verbose:
        print "Loading count table: ", input_count_table
    if (ext == '.gz'):
        genome_table = parse_biom_table(gzip.open(input_count_table,'rb'))
    else:
        genome_table = parse_biom_table(open(input_count_table,'U'))
    if opts.verbose:
        print "Predicting the metagenome..."
    
    partitioned_metagenomes = partition_metagenome_contributions(otu_table,genome_table,limit_to_functions=limit_to_functions)
    output_text = "\n".join(["\t".join(map(str,i)) for i in partitioned_metagenomes])
    if opts.verbose:
        print "Writing results to output file: ",opts.output_fp
        
    make_output_dir_for_file(opts.output_fp)
    open(opts.output_fp,'w').write(output_text)

if __name__ == "__main__":
    main()
