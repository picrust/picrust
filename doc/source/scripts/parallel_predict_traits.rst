.. _parallel_predict_traits:

.. index:: parallel_predict_traits.py

*parallel_predict_traits.py* -- Runs predict_traits.py in parallel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`parallel_predict_traits.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-observed_trait_table
		The input trait table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format
	-t, `-`-tree
		The full reference tree, in Newick format
	-o, `-`-output_trait_table
		The output filepath for trait predictions
	
	**[OPTIONAL]**
		
	-a, `-`-calculate_accuracy_metrics
		If specified, calculate accuracy metrics (i.e. how accurate does PICRUSt expect its predictions to be?) and add to output file [default: False]
	-r, `-`-reconstructed_trait_table
		The input trait table describing reconstructed traits (from `ancestral_state_reconstruction.py <./ancestral_state_reconstruction.html>`_) in tab-delimited format [default: None]
	`-`-output_precalc_file_in_biom
		Instead of outputting the precalculated file in tab-delimited format (with otu ids as row ids and traits as columns) output the data in biom format (with otu as SampleIds and traits as ObservationIds) [default: False]
	-c, `-`-reconstruction_confidence
		The input trait table describing confidence intervals for reconstructed traits (from `ancestral_state_reconstruction.py <./ancestral_state_reconstruction.html>`_) in tab-delimited format [default: None]
	-j, `-`-parallel_method
		Method for parallelizaation. Valid choices are: sge, torque, multithreaded [default: multithreaded]
	-n, `-`-num_jobs
		Number of jobs to be submitted. [default: 2]
	-d, `-`-delay
		Number of seconds to pause between launching each job [default: 0]


**Output:**




Basic

::

	parallel_predict_traits.py -i trait_table.tab -t reference_tree.newick -r asr_counts.tab -c asr_ci.tab -o predict_traits.tab


