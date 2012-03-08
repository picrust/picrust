#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The PICRUST project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from biom.parse import parse_biom_table
from picrust.parse import parse_marker_gene_copy_numbers

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","Normalize the counts in raw_otu_table.biom by dividing by marker gene copy numbers provided in copy_numbers.txt. Write the resulting table to normalized_otu_table.biom.","%prog -i raw_otu_table.biom -c copy_numbers.txt -o normalized_otu_table.biom")]
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--input_otu_fp',type="existing_filepath",help='the input otu table filepath in biom format'),
 make_option('-c','--input_count_fp',type="existing_filepath",help='the input marker gene counts on per otu basis'),
 make_option('-o','--output_otu_fp',type="new_filepath",help='the output otu table filepath in biom format'),
]
script_info['optional_options'] = [
 make_option('--metadata_identifer',
             default='CopyNumber',
             help='identifier for copy number entry as observation metadata [default: %default]'),
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    table = parse_biom_table(open(opts.input_otu_fp,'U'))
    copy_numbers = parse_marker_gene_copy_numbers(open(opts.input_count_fp),
                                                  metadata_identifier=opts.metadata_identifer)
    #Need to only keep data relevant to our otu list
    ids=[]
    for x in table.iterObservations():
        ids.append(str(x[1]))

    copy_numbers_filtered={}
    for x in ids:
        copy_numbers_filtered[x]=copy_numbers[x]

    table.addObservationMetadata(copy_numbers_filtered)

    normalized_table = table.normObservationByMetadata(opts.metadata_identifer)
    open(opts.output_otu_fp,'w').write(\
     normalized_table.getBiomFormatJsonString('PICRUST'))


if __name__ == "__main__":
    main()
