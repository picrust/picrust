.. _run_genome_evaluations:

.. index:: run_genome_evaluations.py

*run_genome_evaluations.py* -- Runs genome evaluations on PICRUSt. 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Using files created by `make_test_datasets.py <./make_test_datasets.html>`_ it runs each test dataset through the ASR (`ancestral_state_reconstruction.py <./ancestral_state_reconstruction.html>`_) and the genome prediction (`predict_traits.py <./predict_traits.html>`_)


**Usage:** :file:`run_genome_evaluations.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_dir
		Directory containing one or more test datasets
	-t, `-`-ref_tree
		Reference tree that was used with make_test_datasets
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		The output directory [default: <input_dir>]
	-j, `-`-parallel_method
		Method for parallelization. Valid choices are: sge, torque, multithreaded [default: multithreaded]
	-m, `-`-prediction_method
		Method for trait prediction.  See `predict_traits.py <./predict_traits.html>`_ for full documentation. Valid choices are: asr_and_weighting, nearest_neighbor, random_neighbor [default: asr_and_weighting]
	`-`-with_confidence
		If set, calculate confidence intervals with ace_ml or ace_reml, and use confidence intervals in trait prediction
	`-`-with_accuracy
		If set, calculate accuracy using the NSTI (nearest sequenced taxon index) during trait prediction
	-a, `-`-asr_method
		Method for ancestral_state_reconstruction.  See `ancestral_state_reconstruction.py <./ancestral_state_reconstruction.html>`_ for full documentation. Valid choices are: ace_ml, ace_reml, ace_pic, wagner [default: wagner]
	-w, `-`-weighting_method
		Method for weighting during trait prediction.  See `predict_traits.py <./predict_traits.html>`_ for full documentation. Valid choices are: linear, exponential, equal [default: exponential]
	-n, `-`-num_jobs
		Number of jobs to be submitted (if --parallel). [default: 100]
	`-`-tmp-dir
		Location to store intermediate files [default: <output_dir>]
	`-`-force
		Run all jobs even if output files exist [default: False]
	`-`-check_for_null_files
		Check if pre-existing output files have null files. If so remove them and re-run. [default: False]


**Output:**

Predictions from `predict_traits.py <./predict_traits.html>`_ for each test dataset.


**Minimum Requirments:**

Provide a directory that contains one or more datasets created by `make_test_datasets.py <./make_test_datasets.html>`_ and the original reference tree used

::

	run_genome_evaluations.py -i test_datasets_dir -t reference_tree_fp

**Specify output file:**

::

	run_genome_evaluations.py -i test_datasets_dir -t reference_tree_fp -o output_dir

**Force the launching of jobs that alredy seem done by overwriting existing output files:**

::

	run_genome_evaluations.py --force -i test_datasets_dir -t reference_tree_fp -o output_dir


