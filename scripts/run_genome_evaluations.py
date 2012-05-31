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
from glob import glob
from picrust.util import make_output_dir_for_file
from cogent.app.util import get_tmp_filename
from picrust.util import get_picrust_project_dir,make_output_dir
from picrust.parallel import submit_jobs, wait_for_output_files
from os.path import join
from re import split

script_info = {}
script_info['brief_description'] = "Runs genome evaluations on PI-CRUST. "
script_info['script_description'] = "\
Using files created by make_test_datasets.py it runs each test dataset through the ASR (ancestral_state_reconstruction.py) and the genome prediction (predict_traits.py)"

script_info['script_usage'] = [\
("Minimum Requirments","Provide a directory that contains one or more datasets created by make_test_datasets.py and the original reference tree used","%prog -i test_datasets_dir -t reference_tree_fp"),\
("Specify output file","","%prog -i test_datasets_dir -t reference_tree_fp -o output_dir")]

script_info['output_description']= "Predictions from predict_traits.py for each test dataset."

script_info['required_options'] = [\
make_option('-i','--input_dir',type="existing_dirpath",help='directory containing one or more test datasets'),\
make_option('-t','--ref_tree',type="existing_filepath",help='reference tree that was used with make_test_datasets'),\
]
parallel_method_choices=['sge','torque','multithreaded']

script_info['optional_options'] = [\
make_option('-o','--output_dir',type="new_dirpath",help='the output directory [default: <input_dir>]'),\
make_option('-j','--parallel_method',type='choice',\
            help='Method for parallelization. Valid choices are: '+\
            ', '.join(parallel_method_choices) + ' [default: %default]',\
            choices=parallel_method_choices,default='multithreaded'),\
make_option('-n','--num_jobs',action='store',type='int',\
            help='Number of jobs to be submitted (if --parallel). [default: %default]',\
            default=100),\
make_option('--tmp-dir',type="new_dirpath",help='location to store intermediate files [default: <output_dir>]'),\

]

script_info['version'] = __version__
       

def main():
    option_parser, opts, args =\
                   parse_command_line_parameters(**script_info)

    #set some defaults for the options
    input_dir=opts.input_dir
    output_dir=opts.output_dir or input_dir
    tmp_dir=opts.tmp_dir or output_dir
    parallel_method=opts.parallel_method


    #create the output directory unless it already exists
    make_output_dir(output_dir)

    if(parallel_method=='sge'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_jobs_sge.py')
    elif(parallel_method=='multithreaded'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_jobs.py')
    elif(parallel_method=='torque'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_jobs_torque.py')
    else:
        raise RuntimeError


    #get the test datasets to run in the input directory (based on exp_traits files)
    expect_test_files=glob(join(input_dir,'exp_traits--*')) 

    test_datasets={}
    for file_name in expect_test_files:
        test_id=file_name.replace(join(input_dir,'exp_traits--'),'',1)
        #create a dict with the test files as values in the ref list
        test_datasets[test_id]=[ join(input_dir,'test_trait_table--'+test_id),join(input_dir,'test_tree--'+test_id),join(input_dir,'exp_traits--'+test_id)]
    
    created_tmp_files=[]    
    output_files=[]

    #create a tmp file to store the job commands (which we will pass to our parallel script to run)
    jobs_fp=get_tmp_filename(tmp_dir=tmp_dir,prefix='jobs_')
    jobs=open(jobs_fp,'w')
    created_tmp_files.append(jobs_fp)

    #get location of scripts we need to run
    asr_script_fp = join(get_picrust_project_dir(),'scripts','ancestral_state_reconstruction.py')
    predict_traits_script_fp = join(get_picrust_project_dir(),'scripts','predict_traits.py')

    asr_method='wagner'

    #run each test dataset through the pipeline
    for test_id in test_datasets:

        asr_out_fp=join(output_dir,'asr--'+test_id)
        created_tmp_files.append(asr_out_fp)
        
        #create the asr command
        asr_cmd= "{0} -i '{1}' -t '{2}' -m {3} -o '{4}'".format(asr_script_fp, test_datasets[test_id][0], test_datasets[test_id][1], asr_method, asr_out_fp)

        predict_traits_out_fp=join(output_dir,'predict_traits--'+test_id)
        output_files.append(predict_traits_out_fp)

        genome_id=split('--',test_id)[2]
        
        #create the predict traits command
        predict_traits_cmd= "{0} -i '{1}' -t '{2}' -r '{3}' -g '{4}' -o '{5}'".format(predict_traits_script_fp, test_datasets[test_id][0], opts.ref_tree, asr_out_fp,genome_id, predict_traits_out_fp)
 
        #add job command to the the jobs file
        jobs.write(asr_cmd+';'+predict_traits_cmd+"\n")

    jobs.close()

    #created_tmp_files.extend(output_files)

    #submit the jobs
    job_prefix='genome_eval_'
    submit_jobs(cluster_jobs_fp ,jobs_fp,job_prefix,num_jobs=opts.num_jobs)

    #wait until all jobs finished (e.g. simple poller)
    wait_for_output_files(output_files)


if __name__ == "__main__":
    main()
