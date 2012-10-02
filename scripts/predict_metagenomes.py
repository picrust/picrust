#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The PICRUST project"
__credits__ = ["Greg Caporaso","Morgan Langille"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from biom.parse import parse_biom_table
from picrust.predict_metagenomes import predict_metagenomes
from picrust.format import format_biom_table
from picrust.util import make_output_dir_for_file
from os import path
import gzip

script_info = {}
script_info['brief_description'] = "This script produces the actual metagenome functional predictions for a given OTU table."
script_info['script_description'] = "This script produces the actual metagenome functional predictions for a given OTU table."
script_info['script_usage'] = [("","Basic usage:","%prog -i otus.biom -c KEGG_acepic_predict_traits_97.biom.gz -o predicted_metagenomes.biom"),
                               ("","Change output format to plain tab-delimited:","%prog -f -i otus.biom -c KEGG_acepic_predict_traits_97.biom.gz -o predicted_metagenomes.tab")]
script_info['output_description']= "Output is a table of function counts (e.g. KEGG KOs) by sample ids."
script_info['required_options'] = [
 make_option('-i','--input_otu_table',type='existing_filepath',help='the input otu table in biom format'),
 make_option('-c','--input_count_table',type="existing_filepath",help='the input trait counts on per otu basis in biom format (can be gzipped)'),
 make_option('-o','--output_metagenome_table',type="new_filepath",help='the output file for the predicted metagenome'),
]
script_info['optional_options'] = [
 make_option('-f','--format_tab_delimited',action="store_true",default=False,help='output the predicted metagenome table in tab-delimited format [default: %default]')
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
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
    predicted_metagenomes = predict_metagenomes(otu_table,genome_table)

    if opts.verbose:
        print "Writing results to output file: ",opts.output_metagenome_table

    make_output_dir_for_file(opts.output_metagenome_table)
    if(opts.format_tab_delimited):
        open(opts.output_metagenome_table,'w').write(predicted_metagenomes.delimitedSelf())
    else:
        open(opts.output_metagenome_table,'w').write(format_biom_table(predicted_metagenomes))

if __name__ == "__main__":
    main()
