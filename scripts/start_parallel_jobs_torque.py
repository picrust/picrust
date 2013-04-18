#!/usr/bin/env python 

"""A simple qsub based cluster submission script for Torque."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project" 
__credits__ = ["Jens Reeder",
               "Rob Knight",
               "Greg Caporaso",
               "Morgan Langille"]
__license__ = "GPL"
__version__ = "0.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from optparse import OptionParser
from os.path import exists
from os import remove, rename, rmdir, makedirs

from cogent.util.misc import app_path
from cogent.app.util import get_tmp_filename
from picrust.make_cluster_jobs import make_torque_jobs, submit_cluster_jobs
from cogent.util.option_parsing import parse_command_line_parameters, make_option

script_info = {}
script_info['brief_description'] = "Starts multiple jobs in parallel on Torque/qsub based multiprocessor systems."
script_info['script_description'] = "This script is designed to start multiple jobs in parallel on cluster systems with a SGE/qsub based scheduling system."
script_info['script_usage'] = [\
 ("Example",\
 "Start each command listed in test_jobs.txt in parallel. The run id for these jobs will be RUNID. ",\
 "%prog -ms test_jobs.txt RUNID")]
script_info['output_description']= "No output is created."
script_info['required_options'] = []
script_info['optional_options'] = [\
    make_option('-m','--make_jobs',action='store_true',\
                    help='make the job files [default: %default]'),

    make_option('-s','--submit_jobs',action='store_true',\
                    help='submit the job files [default: %default]'),

    make_option('-q','--queue',action='store',\
                    type='string',dest='queue', \
                    help='name of queue to submit to '+\
                    ' [default: %default]'),

    make_option('-j','--job_dir', action='store',\
                    type='string',dest='job_dir',\
                    help='directory to store the jobs '+\
                    '[default: %default]', default="jobs/"),
    make_option('-n','--num_jobs',action='store',type='int',\
                help='Number of jobs to group commands into. [default: %default]',\
                default=100)
]

script_info['version'] = __version__
script_info['disallow_positional_arguments'] = False

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
       
    if opts.submit_jobs and not opts.make_jobs:
        option_parser.error('Must pass -m if passing -s. (Sorry about this, '+\
        'it\'s for backwards-compatibility.)') 

    min_args = 2
    if len(args) != min_args:
        option_parser.error('Program requires <commands file> and  <job prefix>')

    if (len(args[1])>10 or len(args[1])==0):
        option_parser.error('job prefix must be 1-10 characters long')
 
    commands = list(open(args[0]))
    job_prefix = args[1]

    if(not exists(opts.job_dir)):
        try:
            makedirs(opts.job_dir)
        except OSError:
            exit(" Jobs directory can not be created. "
                 +"Check for permissions or file with the same name: %s\n"
                 % opts.job_dir)

    if (opts.make_jobs):
        filenames = make_torque_jobs(commands, job_prefix, opts.queue, opts.job_dir,opts.num_jobs)
    else:
        exit("Should we ever get here???")
    if (opts.submit_jobs):
        submit_cluster_jobs(filenames, opts.verbose)

if __name__ == "__main__":
    main()
