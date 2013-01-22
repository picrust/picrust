#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Greg Caporaso","Jesse Zaneveld","Morgan Langille"]
__license__ = "GPL"
__version__ = "0.9.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from cStringIO import StringIO
from cogent.util.option_parsing import parse_command_line_parameters, make_option
from biom.parse import parse_biom_table
from picrust.predict_metagenomes import predict_metagenomes, calc_nsti,\
  load_subset_from_biom_str
from picrust.util import make_output_dir_for_file,format_biom_table
from os import path
from os.path import join
from picrust.util import get_picrust_project_dir
import gzip

script_info = {}
script_info['brief_description'] = "This script produces the actual metagenome functional predictions for a given OTU table."
script_info['script_description'] = ""
script_info['script_usage'] = [("","Predict metagenomes from genomes.biom and otus.biom.","%prog -i normalized_otus.biom -o predicted_metagenomes.biom"),
                               ("","Change output format to plain tab-delimited:","%prog -f -i normalized_otus.biom -o predicted_metagenomes.tab")]
script_info['output_description']= "Output is a table of function counts (e.g. KEGG KOs) by sample ids."
script_info['required_options'] = [
 make_option('-i','--input_otu_table',type='existing_filepath',help='the input otu table in biom format'),
 make_option('-o','--output_metagenome_table',type="new_filepath",help='the output file for the predicted metagenome')
]
script_info['optional_options'] = [\
    make_option('-c','--input_count_table',default=join(get_picrust_project_dir(),'picrust','data','ko_precalculated.biom.gz'),type="existing_filepath",help='Precalculated function predictions on per otu basis in biom format (can be gzipped) [default: %default]'),
    make_option('-a','--accuracy_metrics',default=None,type="new_filepath",help='If provided, calculate accuracy metrics for the predicted metagenome.  NOTE: requires that per-genome accuracy metrics were calculated using predict_traits.py during genome prediction (e.g. there are "NSTI" values in the genome .biom file metadata)'),
    make_option('--suppress_subset_loading',default=False,action="store_true",help='Normally, only counts for OTUs present in the sample are loaded.  If this flag is passed, the full biom table is loaded.  This makes no difference for the analysis, but may result in faster load times (at the cost of more memory usage)'),
  make_option('-f','--format_tab_delimited',action="store_true",default=False,help='output the predicted metagenome table in tab-delimited format [default: %default]')]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    if opts.verbose:
        print "Loading otu table: ",opts.input_otu_table

    otu_table = parse_biom_table(open(opts.input_otu_table,'U'))
    ext=path.splitext(opts.input_count_table)[1]
   
    if not opts.suppress_subset_loading:
        #Now we want to use the OTU table information
        #to load only rows in the count table corresponding
        #to relevant OTUs
        ids_to_load = otu_table.ObservationIds

        if opts.verbose:
            print "Loading traits for %i organisms from the trait table" %len(ids_to_load)
    

    if opts.verbose:
        print "Loading count table: ", opts.input_count_table
    #if (ext == '.gz'):
    #    genome_table = parse_biom_table(gzip.open(opts.input_count_table,'rb'))
    #else:
    #    genome_table = parse_biom_table(open(opts.input_count_table,'U'))

    if (ext == '.gz'):
        genome_table_str = gzip.open(opts.input_count_table,'rb').read()
    else:
        genome_table_str = open(opts.input_count_table,'U').read()
    
    #In the genome/trait table genomes are the samples and 
    #genes are the observations
    
    if not opts.suppress_subset_loading:
        if opts.verbose:
            print "Building a biom Table object from relevant entries"
        genome_table = load_subset_from_biom_str(genome_table_str,ids_to_load,axis='samples')
    else:
        if opts.verbose:
            print "Loading *full* count table because --suppress_subset_loading was passed. This may result in high memory usage"
        genome_table = parse_biom_table(genome_table_str)
    
    if opts.verbose:
        print "Done loading trait table with %i samples" %(len(genome_table.SampleIds))

    make_output_dir_for_file(opts.output_metagenome_table)

    if opts.accuracy_metrics:
        # Calculate accuracy metrics
        #unweighted_nsti = calc_nsti(otu_table,genome_table,weighted=False)
        #print "Unweighted NSTI:", unweighted_nsti
        
        weighted_nsti = calc_nsti(otu_table,genome_table,weighted=True)
        samples= weighted_nsti[0]
        nstis = list(weighted_nsti[1])
        #print "Samples:",samples
        #print "NSTIs:",nstis
        samples_and_nstis = zip(samples,nstis)
        #print "Samples and NSTIs:",samples_and_nstis
        lines = ["#Sample\tMetric\tValue\n"]
        #print weighted_nsti
        for sample,nsti in samples_and_nstis:
            line = "%s\tWeighted NSTI\t%s\n" %(sample,str(nsti))
            lines.append(line)

        if opts.verbose:
            for l in sorted(lines):
                print l
        if opts.verbose:
            print "Writing accuracy information to file:", opts.accuracy_metrics
        open(opts.accuracy_metrics,'w').writelines(sorted(lines))

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
