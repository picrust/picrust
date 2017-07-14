#!/usr/bin/env python
# File created on 1 Feb 2012
from __future__ import division

__author__ = "Morgan Langille"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Morgan Langille","Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.1.1-dev"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from glob import glob
from picrust.util import make_output_dir_for_file,file_contains_nulls
from cogent.app.util import get_tmp_filename
from picrust.util import get_picrust_project_dir,make_output_dir
from picrust.parallel import submit_jobs, wait_for_output_files
from os.path import join,exists
from os import remove
from re import split

script_info = {}
script_info['brief_description'] = "Runs genome evaluations on PICRUSt. "
script_info['script_description'] = "\
Using files created by make_test_datasets.py it runs each test dataset through the ASR (ancestral_state_reconstruction.py) and the genome prediction (predict_traits.py)"

script_info['script_usage'] = [\
("Minimum Requirments","Provide a directory that contains one or more datasets created by make_test_datasets.py and the original reference tree used","%prog -i test_datasets_dir -t reference_tree_fp"),\
("Specify output file","","%prog -i test_datasets_dir -t reference_tree_fp -o output_dir"),\
("Force the launching of jobs that alredy seem done by overwriting existing output files","", "%prog --force -i test_datasets_dir -t reference_tree_fp -o output_dir"),\
]

script_info['output_description']= "Predictions from predict_traits.py for each test dataset."

script_info['required_options'] = [\
make_option('-i','--input_dir',type="existing_dirpath",help='directory containing one or more test datasets'),\
make_option('-t','--ref_tree',type="existing_filepath",help='reference tree that was used with make_test_datasets'),\
]

# Choices for choice options
parallel_method_choices=['sge','torque','multithreaded']
predict_traits_choices =['asr_and_weighting','nearest_neighbor','random_neighbor']
asr_choices = ['ace_ml', 'ace_reml', 'ace_pic', 'wagner']
weighting_choices = ['linear','exponential','equal']

script_info['optional_options'] = [\
make_option('-o','--output_dir',type="new_dirpath",help='the output directory [default: <input_dir>]'),\
make_option('-j','--parallel_method',type='choice',\
            help='Method for parallelization. Valid choices are: '+\
            ', '.join(parallel_method_choices) + ' [default: %default]',\
            choices=parallel_method_choices,default='multithreaded'),\
make_option('-m','--prediction_method',type='choice',\
            help='Method for trait prediction.  See predict_traits.py for full documentation. Valid choices are: '+\
            ', '.join(predict_traits_choices) + ' [default: %default]',\
            choices=predict_traits_choices,default='asr_and_weighting'),\
make_option('--with_confidence',action='store_true',default=False,\
            help='If set, calculate confidence intervals with ace_ml or ace_reml, and use confidence intervals in trait prediction'),\
make_option('--with_accuracy',action='store_true',default=False,\
            help='If set, calculate accuracy using the NSTI (nearest sequenced taxon index) during trait prediction'),\
make_option('-a','--asr_method',type='choice',\
            help='Method for ancestral_state_reconstruction.  See ancestral_state_reconstruction.py for full documentation. Valid choices are: '+\
            ', '.join(asr_choices) + ' [default: %default]',\
            choices=asr_choices,default='wagner'),\
make_option('-w','--weighting_method',type='choice',\
            help='Method for weighting during trait prediction.  See predict_traits.py for full documentation. Valid choices are: '+\
            ', '.join(weighting_choices) + ' [default: %default]',\
            choices=weighting_choices,default='exponential'),\
make_option('-n','--num_jobs',action='store',type='int',\
            help='Number of jobs to be submitted. [default: %default]',\
            default=100),\
make_option('--tmp-dir',type="new_dirpath",help='location to store intermediate files [default: <output_dir>]'),\
make_option('--force',action='store_true',default=False, help='run all jobs even if output files exist [default: %default]'),\
make_option('--check_for_null_files',action='store_true',default=False, help='check if pre-existing output files have null files. If so remove them and re-run. [default: %default]')
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
    asr_method = opts.asr_method
    predict_traits_method = opts.prediction_method
    
    if opts.num_jobs > 20 and parallel_method == 'multithreaded':
        raise ValueError('You probably dont want to run multithreaded evaluations with a large num_jobs. Please adjust options num_jobs and or parallel_method')
        
    if opts.with_confidence and asr_method not in ['ace_ml','ace_reml']:
        raise ValueError("PICRUST currently only supports confidence intervals with the ace_ml and ace_reml ASR methods")

    if opts.verbose:
        print "Reconstruction method:",asr_method
        print "Prediction method:",predict_traits_method
        print "Parallel method:",parallel_method
        print "num_jobs:",opts.num_jobs
        print "\nOutput will be saved here:'%s'" %output_dir 
    
    #create the output directory unless it already exists
    make_output_dir(output_dir)

    if(parallel_method=='sge'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_picrust_jobs_sge.py')
    elif(parallel_method=='multithreaded'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_picrust_jobs.py')
    elif(parallel_method=='torque'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_picrust_jobs_torque.py')
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

    #run each test dataset through the pipeline
    for test_id in test_datasets:

        asr_out_fp=join(output_dir,'asr--'+asr_method+'--'+test_id)
        asr_params_out_fp=join(output_dir,'--'.join(['asr',asr_method,'asr_params',test_id]))
        created_tmp_files.append(asr_out_fp)

        if opts.check_for_null_files and exists(asr_out_fp) and file_contains_nulls(asr_out_fp):
            #remove file
            if opts.verbose:
                print "Existing ASR file contains null characters. Will run ASR again after removing: "+asr_out_fp
            remove(asr_out_fp)
        

        if exists(asr_out_fp) and not opts.force:
            if opts.verbose:
                print "Output file: {0} already exists, so we will skip it.".format(asr_out_fp)
            asr_cmd = "echo 'Skipping ASR for %s, file %s exists already'" %(test_id,asr_out_fp)
        else:
            #create the asr command
            asr_cmd= """python {0} -i "{1}" -t "{2}" -m {3} -o "{4}" -c "{5}" """.format(asr_script_fp, test_datasets[test_id][0], test_datasets[test_id][1], asr_method, asr_out_fp, asr_params_out_fp)

        predict_traits_out_fp=join(output_dir,'--'.join(['predict_traits',predict_traits_method,\
          opts.weighting_method,test_id]))
        
        if opts.with_accuracy:
            predict_traits_accuracy_out_fp=join(output_dir,'--'.join(['predict_traits',predict_traits_method,\
              opts.weighting_method,'accuracy_metrics',test_id]))

        if opts.check_for_null_files and exists(predict_traits_out_fp) and file_contains_nulls(predict_traits_out_fp):
            if opts.verbose:
                print "Existing trait predictions file contains null characters. Will run it again after removing: "+predict_traits_out_fp
            remove(predict_traits_out_fp)

        if exists(predict_traits_out_fp) and not opts.force:
            if opts.verbose:
                print "Prediction file: {0} already exists. Skipping ASR and prediction for this organism".format(predict_traits_out_fp)
            continue
        
        output_files.append(predict_traits_out_fp)

        genome_id=split('--',test_id)[2]
        
        if predict_traits_method == 'nearest_neighbor':
            #don't do asr step
            predict_traits_cmd= """python {0} -i "{1}" -t "{2}" -g "{3}" -o "{4}" -m "{5}" """.format(predict_traits_script_fp, test_datasets[test_id][0], opts.ref_tree, genome_id, predict_traits_out_fp,predict_traits_method)
            jobs.write(predict_traits_cmd+"\n")
        else:

            #create the predict traits command
            predict_traits_cmd= """python {0} -i "{1}" -t "{2}" -r "{3}" -g "{4}" -o "{5}" -m "{6}" -w {7} """.format(predict_traits_script_fp,\
            test_datasets[test_id][0], opts.ref_tree, asr_out_fp,genome_id, predict_traits_out_fp,predict_traits_method,opts.weighting_method)

            #Instruct predict_traits to use confidence intervals output by ASR
            if opts.with_confidence:
                confidence_param = ' -c "%s"' %(asr_params_out_fp)
                predict_traits_cmd = predict_traits_cmd + confidence_param
        
            #Instruct predict traits to output the NTSI measure of distance to
            #nearby sequences.

            if opts.with_accuracy:
                accuracy_param = ' -a "%s"' %(predict_traits_accuracy_out_fp)
                predict_traits_cmd = predict_traits_cmd + accuracy_param

        

 
            #add job command to the the jobs file
            jobs.write(asr_cmd+';'+predict_traits_cmd+"\n")

    jobs.close()

    #created_tmp_files.extend(output_files)

    #submit the jobs
    job_prefix='eval_'
    
    if opts.verbose:
        print "Submitting jobs:",cluster_jobs_fp,jobs_fp,job_prefix,opts.num_jobs
    submit_jobs(cluster_jobs_fp ,jobs_fp,job_prefix,num_jobs=opts.num_jobs)
    #wait until all jobs finished (e.g. simple poller)
    #NOTE: Commented out since nothing happens afterwards currently
    #wait_for_output_files(output_files)


if __name__ == "__main__":
    main()
