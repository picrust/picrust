.. _evaluate_test_datasets:

.. index:: evaluate_test_datasets.py

*evaluate_test_datasets.py* -- Evaluate the accuracy of character predictions, given directories of expected vs. observed test results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The script finds all paired expected and observed values in a set of directories and generates the following output: 1) data for a scatterplot of observed vs. expected values for each character (typically gene family count) within each organism (so one file per organism). 2) A summary of accuracy across all organisms.   
    character 


**Usage:** :file:`evaluate_test_datasets.py [options]`

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
		Pass comma-separated categories to pool results by those metadata categories. Valid categories are: holdout_method, prediction_method,weighting_method,distance and organism. For example, pass "distance" to output results pooled by holdout distance in addition to holdout method and prediction method  [default: distance]


**Output:**

Outputs will be obs,exp data points for the comparison


Evaluate the accuracy of all predictions in a folder, and output summary statistics.

::

	evaluate_test_datasets.py -i obs_otu_table.biom -e exp_otu_table.txt -o./evaluation_results/


