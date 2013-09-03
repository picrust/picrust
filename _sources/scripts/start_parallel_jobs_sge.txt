.. _start_parallel_jobs_sge:

.. index:: start_parallel_jobs_sge.py

*start_parallel_jobs_sge.py* -- Starts multiple jobs in parallel on SGE/qsub based multiprocessor systems.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script is designed to start multiple jobs in parallel on cluster systems with a SGE/qsub based scheduling system.


**Usage:** :file:`start_parallel_jobs_sge.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-m, `-`-make_jobs
		Make the job files [default: None]
	-s, `-`-submit_jobs
		Submit the job files [default: None]
	-d, `-`-delay
		Number of seconds to pause between launching each job [default: 0]
	-q, `-`-queue
		Name of queue to submit to  [default: None]
	-j, `-`-job_dir
		Directory to store the jobs [default: jobs/]
	-n, `-`-num_jobs
		Number of jobs to group commands into. [default: 100]


**Output:**

No output is created.


**Example:**

Start each command listed in test_jobs.txt in parallel. The run id for these jobs will be RUNID. 

::

	start_parallel_jobs_sge.py -ms test_jobs.txt RUNID


