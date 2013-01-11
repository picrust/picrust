#!/usr/bin/env python
# File created on 1 Feb 2012
from __future__ import division

__author__ = "Morgan Langille"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from picrust.count import wagner_for_picrust
from picrust.ace import ace_for_picrust
from picrust.ancestral_state_reconstruction import run_asr_in_parallel
from picrust.util import make_output_dir_for_file,make_output_dir

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
asr_method_choices=['ace_ml','ace_reml','ace_pic','wagner']
parallel_method_choices=['sge','torque','multithreaded']

script_info['optional_options'] = [\
make_option('-m','--asr_method',type='choice',
                help='Method for ancestral state reconstruction. Valid choices are: '+\
                ', '.join(asr_method_choices) + ' [default: %default]',\
                choices=asr_method_choices,default='ace_pic'),\
make_option('-o','--output_fp',type="new_filepath",help='output trait table [default:%default]',default='asr_counts.tab'),\
make_option('-c','--output_ci_fp',type="new_filepath",help='output table containing 95% confidence intervals, loglik, and brownian motion parameters for each asr prediction [default:%default]',default='asr_ci.tab'),\
make_option('-p','--parallel',action="store_true",help='allow parallelization of asr',default=False),\
make_option('-j','--parallel_method',type='choice',
                help='Method for parallelizaation. Valid choices are: '+\
                ', '.join(parallel_method_choices) + ' [default: %default]',\
                choices=parallel_method_choices,default='sge'),\
make_option('-n','--num_jobs',action='store',type='int',\
                help='Number of jobs to be submitted (if --parallel). [default: %default]',\
                default=100),\
make_option('-d','--debug',action="store_true",help='To aid with debugging; get the command that the app controller is going to run',default=False),\
]

script_info['version'] = __version__
       

def main():
    option_parser, opts, args =\
                   parse_command_line_parameters(**script_info)
    
    if(opts.parallel):
        tmp_dir='jobs/'
        make_output_dir(tmp_dir)
        asr_table, ci_table =run_asr_in_parallel(tree=opts.input_tree_fp,table=opts.input_trait_table_fp,asr_method=opts.asr_method,parallel_method=opts.parallel_method, num_jobs=opts.num_jobs,tmp_dir=tmp_dir,verbose=opts.verbose)
    else:
        #call the apporpriate ASR app controller 
        if(opts.asr_method == 'wagner'):
            asr_table = wagner_for_picrust(opts.input_tree_fp,opts.input_trait_table_fp,HALT_EXEC=opts.debug)
        elif(opts.asr_method == 'bayestraits'):
            pass
        elif(opts.asr_method == 'ace_ml'):
            asr_table,ci_table = ace_for_picrust(opts.input_tree_fp,opts.input_trait_table_fp,'ML',HALT_EXEC=opts.debug)
        elif(opts.asr_method == 'ace_pic'):
            asr_table,ci_table = ace_for_picrust(opts.input_tree_fp,opts.input_trait_table_fp,'pic',HALT_EXEC=opts.debug)
        elif(opts.asr_method == 'ace_reml'):
            asr_table,ci_table = ace_for_picrust(opts.input_tree_fp,opts.input_trait_table_fp,'REML',HALT_EXEC=opts.debug)


    #output the table to file
    make_output_dir_for_file(opts.output_fp)
    asr_table.writeToFile(opts.output_fp,sep='\t')

    #output the CI file (unless the method is wagner)
    if not (opts.asr_method == 'wagner'):
        make_output_dir_for_file(opts.output_ci_fp)
        ci_table.writeToFile(opts.output_ci_fp,sep='\t')
        
    

if __name__ == "__main__":
    main()
