#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Greg Caporaso","Jesse Zaneveld","Morgan Langille"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from biom.parse import parse_biom_table
from picrust.predict_metagenomes import predict_metagenomes, calc_nsti
from picrust.util import make_output_dir_for_file,format_biom_table
from os import path
from numpy import around
import gzip

script_info = {}
script_info['brief_description'] = "This script converts metagenomic relative abundance back to sequence counts, by scaling the relative abundnace of each gene in each sample in a biom file by a user-supplied sequencing depth"
script_info['script_description'] = ""
script_info['script_usage'] = [("","Predict metagenomes from genomes.biom and otus.biom.","%prog -i otus.biom -c KEGG_acepic__predict_traits_97.biom.gz -o predicted_metagenomes.biom"),
                               ("","Change output format to plain tab-delimited:","%prog -f -i otus.biom -c KEGG_acepic_predict_traits_97.biom.gz -o predicted_metagenomes.tab")]
script_info['output_description']= "Output is a table of function counts (e.g. KEGG KOs) by sample ids."
script_info['required_options'] = [
 make_option('-s','--input_seq_depth_file',type='existing_filepath',help='an input tab-delimited table, with samples as the first column and an integer sequencing depth as the second'),
 make_option('-i','--input_count_table',type="existing_filepath",help='the input trait counts on per otu basis in biom format (can be gzipped)'),
 make_option('-o','--output_metagenome_table',type="new_filepath",help='the output file for the scaled metagenome')
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    if opts.verbose:
        print "Loading sequencing depth table: ",opts.input_seq_depth_file
    scaling_factors = {}
    for sample_id,depth in parse_seq_count_file(open(opts.input_seq_depth_file,'U')):
        scaling_factors[sample_id]=depth    
    
    ext=path.splitext(opts.input_count_table)[1]

    if opts.verbose:
        print "Loading count table: ", opts.input_count_table
    if (ext == '.gz'):
        genome_table = parse_biom_table(gzip.open(opts.input_count_table,'rb'))
    else:
        genome_table = parse_biom_table(open(opts.input_count_table,'U'))
    
    if opts.verbose:
        print "Scaling the metagenome..."
        
    scaled_metagenomes = scale_metagenomes(genome_table,scaling_factors)

    if opts.verbose:
        print "Writing results to output file: ",opts.output_metagenome_table
        
    make_output_dir_for_file(opts.output_metagenome_table)
    open(opts.output_metagenome_table,'w').write(format_biom_table(scaled_metagenomes))




def make_sample_transformer(scaling_factors):
    def transform_sample(sample_value,sample_id,sample_metadata):
        #print "Sample ",sample_id," value:",sample_value
        scaling_factor = scaling_factors[sample_id]
        #print "Scaling factor:",scaling_factor
        new_val = around(sample_value * scaling_factor)
        #print "New value:",new_val
        return new_val
    return transform_sample

def scale_metagenomes(metagenome_table,scaling_factors):
    """ scale metagenomes from metagenome table and scaling factors 
    """
    #print metagenome_table.SampleIds
    transform_sample_f = make_sample_transformer(scaling_factors)
    new_metagenome_table = metagenome_table.transformSamples(transform_sample_f)
    return new_metagenome_table

def parse_seq_count_file(lines):
    """Extract sample name, counts from seq count file"""
    
    for line in lines:
        sample,depth = line.split("\t")
        depth = int(depth)
        yield sample,depth


if __name__ == "__main__":
    main()
