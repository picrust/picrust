.. _make_test_datasets:

.. index:: make_test_datasets.py

*make_test_datasets.py* -- Generates test datasets for cross-validation studies of PICRUSt's accuracy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`make_test_datasets.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_trait_table
		The input trait table.
	-t, `-`-input_tree
		The input tree in Newick format
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		The output directory.  Duplicate trees, trait tables, expected values and prediction files will be saved here.[default:./test_datasets/]
	`-`-min_dist
		The minimum phylogenetic distance to use with the holdout method, if applicable.  Usually 0.0.[default:0.0]
	`-`-suppress_tree_modification
		If passed, modify only the trait table, not the tree . [default: False]
	`-`-dist_increment
		The phylogenetic distance increment to use with the holdout method, if applicable.[default:0.03]
	`-`-max_dist
		The maximum phylogenetic distance to use with the holdout method, if applicable.[default:0.45]
	`-`-limit_to_tips
		If specified, limit test dataset generation to specified tips (comma-separated).[default:]
	-m, `-`-method
		The test method to use in generating test data.  Valid choices are:exclude_tips_by_distance,randomize_tip_labels_by_distance,collapse_tree_by_distance [default: exclude_tips_by_distance]


**Output:**




Generate holdout test trees from genome_tree.newick, and save results in the directory ./test_holdout_trees/.

::

	make_test_datasets.py -t genome_tree.newick -o ./test_holdout_trees


