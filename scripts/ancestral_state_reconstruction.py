#!/usr/bin/env python
# File created on 1 Feb 2012
from __future__ import division

__author__ = "Morgan Langille"
__copyright__ = "Copyright 2012, The PI-CRUST Project"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from picrust.count import wagner_for_picrust
from picrust.ace import ace_for_picrust

script_info = {}
script_info['brief_description'] = "Runs ancestral state reconstruction given a tree and trait table"
script_info['script_description'] = "\
Provides a common interface for running various ancenstral state reconstruction methods (e.g. ACE, BayesTraits, etc.)."
script_info['script_usage'] = [\
("Minimum Requirments","Provide a tree file and trait table file","%prog -i trait_table_fp -t tree_fp"),\
("Specify output file","","%prog -i trait_table_fp -t tree_fp -o output_file_fp")]
script_info['output_description']= "A table containing trait information for internal nodes of the tree."

script_info['required_options'] = [\
make_option('-t','--input_tree_fp',type="existing_filepath",help='the tree to use for ASR'),\
make_option('-i','--input_trait_table_fp',type="existing_filepath",help='the trait table to use for ASR'),\
]
asr_method_choices=['bayestraits','ace_ml','ace_reml','ace_pic','wagner']
script_info['optional_options'] = [\
make_option('-m','--asr_method',type='choice',
                help='Method for ancestral state reconstruction. Valid choices are: '+\
                ', '.join(asr_method_choices) + ' [default: %default]',\
                choices=asr_method_choices,default='wagner'),\
make_option('-o','--output_fp',type="new_filepath",help='the output trait table'),\
make_option('-p','--output_prob_fp',type="new_filepath",help='the probablities for each ASR [default:%default]',default='output_file_prob'),\
]

script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    #call the apporpriate ASR app controller 
    if(opts.asr_method == 'wagner'):
        asr_table = wagner_for_picrust(opts.input_tree_fp,opts.input_trait_table_fp)
    elif(opts.asr_method == 'bayestraits'):
        pass
    elif(opts.asr_method == 'ace_ml'):
        asr_table = ace_for_picrust(opts.input_tree_fp,opts.input_trait_table_fp,'ML')
    elif(opts.asr_method == 'ace_pic'):
        asr_table = ace_for_picrust(opts.input_tree_fp,opts.input_trait_table_fp,'pic')
    elif(opts.asr_method == 'ace_reml'):
        asr_table = ace_for_picrust(opts.input_tree_fp,opts.input_trait_table_fp,'REML')
    
    if opts.output_fp:
        asr_table.writeToFile(opts.output_fp,sep='\t')
    else:
        print asr_table.tostring(sep='\t')


if __name__ == "__main__":
    main()
