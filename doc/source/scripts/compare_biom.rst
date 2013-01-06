.. _compare_biom:

.. index:: compare_biom.py

*compare_biom.py* -- Compare the accuracy of biom files (expected and observed) either by observations (default) or by samples.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

 


**Usage:** :file:`compare_biom.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-e, `-`-exp_trait_table_fp
		The expected trait table (biom format)
	-o, `-`-output_fp
		The output file
	
	**[OPTIONAL]**
		
	-c, `-`-compare_observations
		Calculate accuracy values by comparing between observations (instead of between samples) [default: False]
	-n, `-`-normalize
		Convert both expected and observed tables to relative abundances (instead of observations) [default: False]
	-l, `-`-limit_to_expected_observations
		Ignore observations that are not in the expected table[default: False]
	`-`-limit_to_observed_observations
		Ignore observations that are not in the observed table[default: False]
	-s, `-`-shuffle_samples
		Shuffle samples ids randomly before measuring accuracy[default: False]
	`-`-not_relative_abundance_scores
		Round numbers (instead of taking ceil() which is used for RA) before calculating TP,FP,FN,TN [default: False]


**Output:**

Outputs will be tab delimited file with various accuracy metrics.


::

	compare_biom.py -e exp_otu_table.biom -o results.tab obs_otu_table.biom [obs_otu_table2.biom]


