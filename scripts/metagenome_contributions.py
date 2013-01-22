#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011, The PICRUST project"
__credits__ = ["Greg Caporaso","Jesse Zaneveld","Morgan Langille"]
__license__ = "GPL"
__version__ = "0.9.0"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"
 


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from biom.parse import parse_biom_table
from picrust.predict_metagenomes import predict_metagenomes, calc_nsti
from picrust.metagenome_contributions import partition_metagenome_contributions
from picrust.util import make_output_dir_for_file
from os import path
import gzip

script_info = {}
script_info['brief_description'] = "This script partitions metagenome functional contributions according to function, OTU, and sample, for a given OTU table."
script_info['script_description'] = ""
script_info['script_usage'] = [("","Partition the predicted contribution to the  metagenomes from each organism in otus.biom, using the predicted genes for each organism in genes.biom.","%prog -i otus.biom -c KEGG_acepic__predict_traits_97.biom.gz -o predicted_metagenomes.biom"),
                               ("","Change output format to plain tab-delimited:","%prog -f -i otus.biom -c KEGG_acepic_predict_traits_97.biom.gz -o predicted_metagenomes.tab")]
script_info['output_description']= "Output is a table of function counts (e.g. KEGG KOs) by sample ids."
script_info['required_options'] = [
 make_option('-i','--input_otu_table',type='existing_filepath',help='the input otu table in biom format'),
 make_option('-c','--input_count_table',type="existing_filepath",help='the input trait counts on per otu basis in biom format (can be gzipped)'),
 make_option('-o','--output_metagenome_table',type="new_filepath",help='the output file for the predicted metagenome')
]
script_info['optional_options'] = [\
        make_option('-a','--accuracy_metrics',default=None,type="new_filepath",help='If provided, calculate accuracy metrics for the predicted metagenome.  NOTE: requires that per-genome accuracy metrics were calculated using predict_traits.py during genome prediction (e.g. there are "NSTI" values in the genome .biom file metadata)'), 
        
        make_option('--limit_to_function',default=None,help='If provided, only output predictions for the specified function ids.  Multiple function ids can be passed using comma delimiters.'),
        make_option('-f','--format_tab_delimited',action="store_true",default=False,help='output the predicted metagenome table in tab-delimited format [default: %default]')]
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
    ext=path.splitext(opts.input_count_table)[1]

    if opts.verbose:
        print "Loading count table: ", opts.input_count_table
    if (ext == '.gz'):
        genome_table = parse_biom_table(gzip.open(opts.input_count_table,'rb'))
    else:
        genome_table = parse_biom_table(open(opts.input_count_table,'U'))
    if opts.verbose:
        print "Predicting the metagenome..."
    
    partitioned_metagenomes = partition_metagenome_contributions(otu_table,genome_table,limit_to_functions=limit_to_functions)
    output_text = "\n".join(["\t".join(map(str,i)) for i in partitioned_metagenomes])
    if opts.verbose:
        print "Writing results to output file: ",opts.output_metagenome_table
        
    make_output_dir_for_file(opts.output_metagenome_table)
    open(opts.output_metagenome_table,'w').write(output_text)

if __name__ == "__main__":
    main()
