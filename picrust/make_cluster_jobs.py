#!/usr/bin/env python 

"""A simple qsub based cluster submission script."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project" 
__credits__ = ["Jens Reeder", "Rob Knight", "Morgan Langille"]#remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.9.2-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Development"

from os.path import exists
from os import remove, rename, rmdir, makedirs
from subprocess import Popen, PIPE, STDOUT

from cogent.util.misc import app_path, create_dir
from cogent.app.util import ApplicationNotFoundError
from cogent.app.util import get_tmp_filename
from picrust.parallel import grouper
from math import ceil
from time import sleep

# qsub template
SGE_QSUB_TEXT = """#!/bin/bash
##
#$ -cwd

# Wall time
#$ -l h_rt=%s

#$ -S /bin/bash
#$ -e %s
#$ -o %s

hostname

%s
"""

#requires format string (walltime, ncpus, nodes, queue, job_name, keep_output, command)
QSUB_TEXT = """# Walltime Limit: hh:nn:ss 
#PBS -l walltime=%s 

# Node Specification:
#PBS -l ncpus=%d -l nodes=%d

# Queue: Defaults to friendlyq 
#PBS -q %s 

# Mail: options are (a) aborted, (b) begins execution, (e) ends execution
# use -M <email> for additional recipients
# supress email notification
#PBS -m n

# Job Name:
#PBS -N %s 

# Keep output
#PBS -k %s

echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------
cd $PBS_O_WORKDIR 
%s
"""


def make_sge_jobs(commands, job_prefix, queue, jobs_dir="jobs/",num_jobs=100,max_hours_per_job=24):
    """prepare qsub text files.
    
    command: list of commands
    
    job_prefix: a short, descriptive name for the job.

    queue: name of the queue to submit to
    
    jobs_dir: path to directory where job submision scripts are written

    max_hours_per_job: the maximum expected time for each command (this will be multiplied by number of commands per job to get a 'walltime' 
    
    ncpus: number of cpus
    
    nodes: number of nodes
    
    keep_output: keep standard error, standard out, both, or neither
                 o=std out, e=std err, oe=both, n=neither
    """

    filenames=[]
    create_dir(jobs_dir)

    #calculate the number of commands to put in each job
    num_commands_per_job=int(ceil(len(commands)/float(num_jobs)))

    #calculate the walltime (time before job will be killed by scheduler if still running)
    total_time = max_hours_per_job*num_commands_per_job
    walltime= "{0}:00:00".format(total_time)
    
    for command_group in grouper(commands,num_commands_per_job,''):
        job_name = get_tmp_filename(tmp_dir=jobs_dir, prefix=job_prefix+"_",
                                    suffix = ".txt")
        out_fh = open(job_name,"w")

        stderr_fp = job_name+"_stderr"
        stdout_fp = job_name+"_stdout"
        out_fh.write(SGE_QSUB_TEXT % (walltime, stderr_fp, stdout_fp, "\n".join(command_group)))        
        out_fh.close()
        filenames.append(job_name)
    return filenames

def make_torque_jobs(commands, job_prefix, queue, jobs_dir="jobs/",
              walltime="72:00:00", ncpus=1, nodes=1, keep_output="oe"):
    """prepare qsub text files.
    
    command: list of commands
    
    job_prefix: a short, descriptive name for the job.

    queue: name of the queue to submit to
    
    jobs_dir: path to directory where job submision scripts are written

    walltime: the maximal walltime 
    
    ncpus: number of cpus
    
    nodes: number of nodes
    
    keep_output: keep standard error, standard out, both, or neither
                 o=std out, e=std err, oe=both, n=neither
    """

    filenames=[]
    create_dir(jobs_dir)
    for command in commands:
        job_name = get_tmp_filename(tmp_dir=jobs_dir, prefix=job_prefix+"_",
                                    suffix = ".txt")
        out_fh = open(job_name,"w")

        out_fh.write(QSUB_TEXT % (walltime, ncpus, nodes, queue, job_prefix,
                                  keep_output, command))        
        out_fh.close()
        filenames.append(job_name)
    return filenames

def submit_cluster_jobs(filenames, verbose=False, delay=0):
    """Submit jobs in filenames.

    filenames: list of prepared qsub job scripts, ready to be submitted

    verbose: a binary verbose flag
    """
    if(not app_path("qsub")):
        raise ApplicationNotFoundError,"qsub not found. Can't submit jobs."
    
    for file in filenames:        
        command = 'qsub %s' % file
        result = Popen(command, shell=True, universal_newlines=True,\
                           stdout=PIPE, stderr=STDOUT).stdout.read()
        if verbose:
            print result
        sleep(delay)
