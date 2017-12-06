.. _pool_test_datasets:

.. index:: pool_test_datasets.py

*pool_test_datasets.py* -- Pool character predictions within a directory, given directories of expected vs. observed test results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The script finds all paired expected and observed values in a set of directories and generates pooled .biom files in a specified output directory


**Usage:** :file:`pool_test_datasets.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-trait_table_dir
		The input trait table directory (files in biom format)
	-e, `-`-exp_trait_table_dir
		The input expected trait table directory (files in biom format)
	-o, `-`-output_dir
		The output directory
	
	**[OPTIONAL]**
		
	-f, `-`-field_order
		Pass comma-separated categories, in the order they appear in file names.   Categories are "file_type","prediction_method","weighting_method","holdout_method" (randomization vs. holdout),"distance",and "organism".  Example:  "-f file_type,test_method,asr_method specifies that files will be in the form: predict_traits--distance_exclusion--wagner.  Any unspecified values are set to "not_specified".  [default: file_type,prediction_method,weighting_method,holdout_method,distance,organism]
	-p, `-`-pool_by
		Pass comma-separated categories to pool results by those metadata categories. Valid categories are: holdout_method, prediction_method,weighting_method,distance and organism. For example, pass "distance" to output results pooled by holdout distance in addition to holdout method and prediction method  [default: False]


**Output:**

Outputs will be obs,exp data points for the comparison


Pool .biom files according to holdout_distance.

::

	pool_test_datasets.py -i obs_otu_table_dir -e exp_otu_table_dir -p distance -o./evaluation_results/pooled_by_distance/


