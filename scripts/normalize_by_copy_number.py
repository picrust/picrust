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
from biom.parse import parse_biom_table, parse_classic_table_to_rich_table
from biom.table import table_factory,DenseOTUTable
from picrust.util import make_output_dir_for_file
from os import path
import gzip

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","Normalize the counts in raw_otu_table.biom by dividing by marker gene copy numbers provided in copy_numbers.biom. Write the resulting table to normalized_otu_table.biom.","%prog -i raw_otu_table.biom -c copy_numbers.biom -o normalized_otu_table.biom")]
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--input_otu_fp',type="existing_filepath",help='the input otu table filepath in biom format'),
 make_option('-c','--input_count_fp',type="existing_filepath",help='the input marker gene counts on per otu basis in biom format'),
 make_option('-o','--output_otu_fp',type="new_filepath",help='the output otu table filepath in biom format'),
]
script_info['optional_options'] = [
 make_option('--metadata_identifer',
             default='CopyNumber',
             help='identifier for copy number entry as observation metadata [default: %default]'),
 make_option('-f','--input_format_classic', action="store_true", default=False, help='input otu table (--input_otu_fp) is in classic Qiime format [default: %default]'),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    if opts.input_format_classic:
        otu_table=parse_classic_table_to_rich_table(open(opts.input_otu_fp,'U'),None,None,None,DenseOTUTable)
    else:
        otu_table = parse_biom_table(open(opts.input_otu_fp,'U'))

    ext=path.splitext(opts.input_count_fp)[1]
    if (ext == '.gz'):
        count_table = parse_biom_table(gzip.open(opts.input_count_fp,'rb'))
    else:
        count_table = parse_biom_table(open(opts.input_count_fp,'U'))
        
    #Need to only keep data relevant to our otu list
    ids=[]
    for x in otu_table.iterObservations():
        ids.append(str(x[1]))

    ob_id=count_table.ObservationIds[0]

    filtered_otus=[]
    filtered_values=[]
    for x in ids:
        if count_table.sampleExists(x):
            filtered_otus.append(x)
            filtered_values.append(otu_table.observationData(x))

    #filtered_values = map(list,zip(*filtered_values))
    filtered_otu_table=table_factory(filtered_values,otu_table.SampleIds,filtered_otus, constructor=DenseOTUTable)

    copy_numbers_filtered={}
    for x in filtered_otus:
        value = count_table.getValueByIds(ob_id,x)
        try:
            #data can be floats so round them and make them integers
            value = int(round(float(value)))
            
        except ValueError:
            raise ValueError,\
                  "Invalid type passed as copy number for OTU ID %s. Must be int-able." % (value)
        if value < 1:
            raise ValueError, "Copy numbers must be greater than or equal to 1."

        copy_numbers_filtered[x]={opts.metadata_identifer:value}
        
    filtered_otu_table.addObservationMetadata(copy_numbers_filtered)
            

    normalized_table = filtered_otu_table.normObservationByMetadata(opts.metadata_identifer)

    make_output_dir_for_file(opts.output_otu_fp)
    open(opts.output_otu_fp,'w').write(\
     normalized_table.getBiomFormatJsonString('PICRUST'))


if __name__ == "__main__":
    main()
