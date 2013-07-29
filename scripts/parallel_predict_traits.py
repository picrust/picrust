#!/usr/bin/env python
from __future__ import division

__author__ = "Morgan Langille"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "0.9.2-dev"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from picrust.util import make_output_dir_for_file,make_output_dir, get_picrust_project_dir, convert_precalc_to_biom
from picrust.format_tree_and_trait_table import load_picrust_tree
from picrust.parallel import submit_jobs, system_call,wait_for_output_files,grouper
from os import makedirs, remove, popen
from os.path import join,splitext
from cogent.app.util import get_tmp_filename

script_info = {}
script_info['brief_description'] = "Runs predict_traits.py in parallel"
script_info['script_description'] = ""
script_info['script_usage'] = [\
("","Basic","%prog -i trait_table.tab -t reference_tree.newick -r asr_counts.tab -c asr_ci.tab -o predict_traits.tab")]
script_info['output_description']= ""

script_info['required_options'] = [\
    make_option('-i','--observed_trait_table',type="existing_filepath",\
                    help='the input trait table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format'),
    make_option('-t','--tree',type="existing_filepath",\
                    help='the full reference tree, in Newick format'),
    make_option('-o','--output_trait_table',type="new_filepath",\
                    help='the output filepath for trait predictions'),\
]

parallel_method_choices=['sge','torque','multithreaded']

script_info['optional_options'] = [\
    make_option('-a','--calculate_accuracy_metrics',default=False,action="store_true",\
                    help='if specified, calculate accuracy metrics (i.e. how accurate does PICRUSt expect its predictions to be?) and add to output file [default: %default]'),
    make_option('-r','--reconstructed_trait_table',
                type="existing_filepath",default=None,
                help='the input trait table describing reconstructed traits (from ancestral_state_reconstruction.py) in tab-delimited format [default: %default]'),

    make_option('--output_precalc_file_in_biom',default=False,action="store_true",
                help='Instead of outputting the precalculated file in tab-delimited format (with otu ids as row ids and traits as columns) output the data in biom format (with otu as SampleIds and traits as ObservationIds) [default: %default]'),
    make_option('-c','--reconstruction_confidence',
                type="existing_filepath",default=None,
                help='the input trait table describing confidence intervals for reconstructed traits (from ancestral_state_reconstruction.py) in tab-delimited format [default: %default]'),

    make_option('-j','--parallel_method',type='choice',
                help='Method for parallelizaation. Valid choices are: '+\
                    ', '.join(parallel_method_choices) + ' [default: %default]',\
                    choices=parallel_method_choices,default='multithreaded'),
    make_option('-n','--num_jobs',action='store',type='int', help='Number of jobs to be submitted. [default: %default]',
                default=2),
    make_option('-d','--delay',action='store',type='int',default=0,
                    help='Number of seconds to pause between launching each job [default: %default]'),

    make_option('--already_calculated',type='existing_filepath',default=None,
                    help='Precalculated file that is missing some otu predictions. Output will contain predictions from this file and the new predictions as well. [default: %default]')
]


script_info['version'] = __version__
       

def combine_predict_trait_output(files):
    #add header
    combined=open(files[0]).readline()
    
    for file_name in files:
        fh=open(file_name)
        #throw away header
        fh.readline()
        combined+=fh.read()
        combined+="\n"

    return combined

def get_tips_not_in_precalc(ids,precalc):
    ids_in_precalc=[]
    for line in open(precalc):
        ids_in_precalc.append(line.strip().split('\t')[0])

    return list(set(ids) - set(ids_in_precalc))


def main():
    option_parser, opts, args =\
                   parse_command_line_parameters(**script_info)

    tmp_dir='jobs/'
    make_output_dir(tmp_dir)

    #Run the jobs
    script_fp = join(get_picrust_project_dir(),'scripts','predict_traits.py')

    if(opts.parallel_method=='sge'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_jobs_sge.py')
    elif(opts.parallel_method=='multithreaded'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_jobs.py')
    elif(opts.parallel_method=='torque'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_jobs_torque.py')
    else:
        raise RuntimeError

    if(opts.verbose):
        print "Loading tree..."
        
    tree = load_picrust_tree(opts.tree, opts.verbose)

    all_tips = [tip.Name for tip in tree.tips()]
    
    if(opts.verbose):
        print "Total number of possible tips to predict: {0}".format(len(all_tips))

    created_tmp_files=[]
    output_files={}
    output_files['counts']=[]
    if opts.reconstruction_confidence:
        output_files['variances']=[]
        output_files['upper_CI']=[]
        output_files['lower_CI']=[]

    if opts.already_calculated:
        all_tips=get_tips_not_in_precalc(all_tips,opts.already_calculated)
        if opts.verbose:
            print "After taking into account tips already predicted, the number of tips left to predict is: {0}".format(len(all_tips))

    #create a tmp file to store the job commands (which we will pass to our parallel script to run)
    jobs_fp=get_tmp_filename(tmp_dir=tmp_dir,prefix='jobs_')
    jobs=open(jobs_fp,'w')
    created_tmp_files.append(jobs_fp)

    if(opts.verbose):
        print "Creating temporary input files in: ",tmp_dir
    
    num_tips_per_job=1000
    for tips_to_predict in [all_tips[i:i+num_tips_per_job] for i in range(0, len(all_tips), num_tips_per_job)]:
        
        #create tmp output files
        tmp_output_fp=get_tmp_filename(tmp_dir=tmp_dir,prefix='out_predict_traits_')
        output_files['counts'].append(tmp_output_fp)

        tip_to_predict_str=','.join(list(tips_to_predict))

        if opts.reconstruction_confidence:
            outfile_base,extension = splitext(tmp_output_fp)
            output_files['variances'].append(outfile_base+"_variances.tab")
            output_files['upper_CI'].append(outfile_base+"_upper_CI.tab")
            output_files['lower_CI'].append(outfile_base+"_lower_CI.tab")
            
            #create the job command
            cmd= "{0} -i {1} -t {2} -r {3} -c {4} -g {5} -o {6}".format(script_fp, opts.observed_trait_table, opts.tree, opts.reconstructed_trait_table, opts.reconstruction_confidence, tip_to_predict_str, tmp_output_fp)

        else:
            cmd= "{0} -i {1} -t {2} -r {3} -g {4} -o {5}".format(script_fp, opts.observed_trait_table, opts.tree, opts.reconstructed_trait_table, tip_to_predict_str, tmp_output_fp)
            

        #NOTE: Calculating NSTI this way is convenient, 
        #but would probably be faster if we ran the NSTI calculation separate (using the --output_accuracy_metrics_only) and added it to the output file later on.
        if opts.calculate_accuracy_metrics:
            cmd=cmd+" -a"

        #add job command to the the jobs file
        jobs.write(cmd+"\n")

    jobs.close()

    #add all output files to tmp list (used later for deletion)
    for predict_type in output_files:
        created_tmp_files.extend(output_files[predict_type])
    if(opts.verbose):
        print "Launching parallel jobs."

    if opts.already_calculated:
        output_files['counts'].append(opts.already_calculated)
        
    #run the job command
    job_prefix='picrust'
    submit_jobs(cluster_jobs_fp ,jobs_fp,job_prefix,num_jobs=opts.num_jobs,delay=opts.delay)

    if(opts.verbose):
        print "Jobs are now running. Will wait until finished."

    #wait until all jobs finished (e.g. simple poller)
    wait_for_output_files(output_files['counts'])

    if(opts.verbose):
        print "Jobs are done running."

    make_output_dir_for_file(opts.output_trait_table)
    outfile_base,extension = splitext(opts.output_trait_table)
    for predict_type in sorted(output_files):
       #Combine output files
        if opts.verbose:
            print "Combining all output files for "+ predict_type

        combined_predictions=combine_predict_trait_output(output_files[predict_type])
        
        if opts.verbose:
            print "Writing combined file for "+predict_type

        if predict_type == 'counts':
        #Output in whatever format the user wants
            if opts.output_precalc_file_in_biom:
                open(opts.output_trait_table,'w').write(format_biom_table(convert_precalc_to_biom(combined_predictions)))
            else:
                open(opts.output_trait_table,'w').write(combined_predictions)
        else:
            if opts.output_precalc_file_in_biom:
                open(outfile_base+"_"+predict_type+".biom",'w').write(format_biom_table(convert_precalc_to_biom(combined_predictions)))
            else:
                open(outfile_base+"_"+predict_type+".tab",'w').write(combined_predictions)    
        
    #clean up all tmp files
    for file in created_tmp_files:
        remove(file)


if __name__ == "__main__":
    main()
