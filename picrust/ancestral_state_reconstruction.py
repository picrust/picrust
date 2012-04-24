#!/usr/bin/env python

from __future__ import division

__author__ = "Morgan Langille"
__copyright__ = "Copyright 2012, The PI-CRUST Project"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"


from os.path import dirname,isdir, exists
from os import makedirs, remove, popen
from subprocess import Popen, PIPE, STDOUT
from picrust.count import wagner_for_picrust
from picrust.ace import ace_for_picrust
from cogent import LoadTable
from cogent.util.table import Table
from cogent.app.util import get_tmp_filename
from picrust.util import get_picrust_project_dir
from os.path import join
from time import sleep



def submit_jobs(path_to_cluster_jobs, jobs_fp, job_prefix,num_jobs=100):
    """ Submit the jobs to the queue using cluster_jobs.py
    """
    cmd = '%s -n %s -ms %s %s' % (path_to_cluster_jobs, num_jobs, jobs_fp, job_prefix)
    stdout, stderr, return_value = system_call(cmd)
    if return_value != 0:
        msg = "\n\n*** Could not start parallel jobs. \n" +\
         "Command run was:\n %s\n" % cmd +\
         "Command returned exit status: %d\n" % return_value +\
         "Stdout:\n%s\nStderr\n%s\n" % (stdout,stderr)
        raise RuntimeError, msg
    
    # Leave this comments in as they're useful for debugging.
    # print 'Return value: %d\n' % return_value
    # print 'STDOUT: %s\n' % stdout
    # print 'STDERR: %s\n' % stderr


def system_call(cmd):
    """ Call cmd and return (stdout, stderr, return_value)"""
    proc = Popen(cmd,shell=True,universal_newlines=True,\
                 stdout=PIPE,stderr=PIPE)
    # communicate pulls all stdout/stderr from the PIPEs to 
    # avoid blocking -- don't remove this line!
    stdout, stderr = proc.communicate()
    return_value = proc.returncode
    return stdout, stderr, return_value
    

def wait_for_output_files(files):
    ''' Function waits until all files exist in the filesystem'''
    #make a copy of the list
    waiting_files=list(files)

    #wait until nothing left in the list
    while(waiting_files):
        #wait 30 seconds between each check
        sleep(30)
        #check each file and keep ones that don't yet exist
        waiting_files=filter(lambda x: not exists(x),waiting_files)


def combine_asr_tables(output_files):
    """ Combine all tables coming from asr output. Cuts 2nd column out and joins them together into single table.
    Assumes all output files have same row identifiers and that these are in the same order.
    """

    #Going to store an array of arrays here
    combined_table=[]
   
    #load in the first column (containing row ids). File doesn't matter since they should all have identical first columns.
    table=LoadTable(filename=output_files[0],header=True,sep='\t')
    row_ids = table.getRawData(columns=[table.Header[0]])
    combined_table.append([table.Header[0]])
    for row_id in row_ids:
        combined_table.append([row_id])
            
    #Now add the rest of the files to the table
    for output_file in output_files:
        #pull out the second column (first column with actual preditions)
        table=LoadTable(filename=output_file,header=True,sep='\t')
        predictions = table.getRawData(columns=[table.Header[1]])

        #Add the header for our column to the list of headers
        combined_table[0].append(table.Header[1])

        #Add rest of values in the column 
        j=1
        for prediction in predictions:
            combined_table[j].append(prediction)
            j+=1

    return combined_table

def run_asr_in_parallel(tree, table, asr_method, parallel_method='sge',tmp_dir='jobs/',num_jobs=100):
    '''Runs the ancestral state reconstructions in parallel'''

    asr_script_fp = join(get_picrust_project_dir(),'scripts','ancestral_state_reconstruction.py')

    if(parallel_method=='sge'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_jobs_sge.py')
    elif(parallel_method=='multithreaded'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_jobs.py')
    elif(parallel_method=='torque'):
        cluster_jobs_fp=join(get_picrust_project_dir(),'scripts','start_parallel_jobs_torque.py')
    else:
        raise RuntimeError
    
    #foreach trait in the table, create a new tmp file with just that trait, and create the job command and add it a tmp jobs file
    table=LoadTable(filename=table, header=True, sep='\t')

    #get dimensions of the table
    dim=table.Shape
    
    created_tmp_files=[]    
    output_files=[]

    #create a tmp file to store the job commands (which we will pass to our parallel script to run)
    jobs_fp=get_tmp_filename(tmp_dir=tmp_dir,prefix='jobs_asr_')
    jobs=open(jobs_fp,'w')
    created_tmp_files.append(jobs_fp)

    #iterate over each column
    for i in range(1,dim[1]):
        #create a new table with only a single trait
        single_col_table=table.getColumns([0,i])
        
        #write the new table to a tmp file
        single_col_fp=get_tmp_filename(tmp_dir=tmp_dir,prefix='in_asr_')
        single_col_table.writeToFile(single_col_fp,sep='\t')
        created_tmp_files.append(single_col_fp)

        #create the job command
        tmp_output_fp=get_tmp_filename(tmp_dir=tmp_dir,prefix='out_asr_')
        output_files.append(tmp_output_fp)
        #tmp_prob_fp=get_tmp_filename()

        cmd= "{0} -i {1} -t {2} -m {3} -o {4}".format(asr_script_fp, single_col_fp, tree, asr_method, tmp_output_fp)
 
        #add job command to the the jobs file
        jobs.write(cmd+"\n")

    jobs.close()
    created_tmp_files.extend(output_files)

    #run the job command
    job_prefix='asr'
    submit_jobs(cluster_jobs_fp ,jobs_fp,job_prefix,num_jobs=num_jobs)

    #wait until all jobs finished (e.g. simple poller)
    wait_for_output_files(output_files)

    #Combine output files
    combined_table=combine_asr_tables(output_files)
    
    #create a Table object
    combined_table=Table(header=combined_table[0],rows=combined_table[1:])
        
    #clean up all tmp files
    for file in created_tmp_files:
        remove(file)

    #return the combined table
    return combined_table
