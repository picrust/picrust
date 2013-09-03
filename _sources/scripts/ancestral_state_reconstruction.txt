.. _ancestral_state_reconstruction:

.. index:: ancestral_state_reconstruction.py

*ancestral_state_reconstruction.py* -- Runs ancestral state reconstruction given a tree and trait table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Provides a common interface for running various ancenstral state reconstruction methods (e.g. ACE, BayesTraits, etc.).


**Usage:** :file:`ancestral_state_reconstruction.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-t, `-`-input_tree_fp
		The tree to use for ASR
	-i, `-`-input_trait_table_fp
		The trait table to use for ASR
	
	**[OPTIONAL]**
		
	-m, `-`-asr_method
		Method for ancestral state reconstruction. Valid choices are: ace_ml, ace_reml, ace_pic, wagner [default: ace_pic]
	-o, `-`-output_fp
		Output trait table [default:asr_counts.tab]
	-c, `-`-output_ci_fp
		Output table containing 95% confidence intervals, loglik, and brownian motion parameters for each asr prediction [default:asr_ci.tab]
	-p, `-`-parallel
		Allow parallelization of asr
	-j, `-`-parallel_method
		Method for parallelizaation. Valid choices are: sge, torque, multithreaded [default: sge]
	-n, `-`-num_jobs
		Number of jobs to be submitted (if --parallel). [default: 100]
	-d, `-`-debug
		To aid with debugging; get the command that the app controller is going to run


**Output:**

A table containing trait information for internal nodes of the tree.


**Example 1:**

Provide a tree file and trait table file:

::

	ancestral_state_reconstruction.py -i trait_table.tab -t pruned_tree.newick -o asr_counts.tab -c asr_ci.tab


