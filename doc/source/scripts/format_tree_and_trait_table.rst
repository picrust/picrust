.. _format_tree_and_trait_table:

.. index:: format_tree_and_trait_table.py

*format_tree_and_trait_table.py* -- Formatting script for filtering and reformatting trees and trait tables.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Reformats scripts and trait tables.  Optional fixes include:  
        -- Add short (epsilon) branch lengths in place of 0 length branches 
        -- Filter out taxa that don't match between tree and trait table 
        -- Output tree in NEXUS format 
        -- Ensure tree is bifurcating (remove polytomies using very short branches)
        -- Convert floating point trait values to integers
        -- Add a short branch length to the root branch (required by BayesTraits)
        -- Remove internal node names (required by BayesTraits)
        


**Usage:** :file:`format_tree_and_trait_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-t, `-`-input_tree
		The input tree (Newick format)
	-i, `-`-input_trait_table
		The input trait table (QIIME OTU table format)
	
	**[OPTIONAL]**
		
	-m, `-`-tree_to_trait_mapping
		A two-column, tab-delimited text file mapping identifiers in the tree(column 1) to identifiers in the trait table (column 2). If supplied, the identifiers in the trait table will be converted to match the identifiers in the tree. (This mapping does not need to be supplied if the tree and trait table already use a common set of identifiers.) [default: None]
	-o, `-`-output_dir
		The output directory [default: ./formatted/]
	`-`-input_table_delimiter
		The character delimiting fields in the input trait table. Valid choices are:tab,space,comma [default: tab]
	`-`-output_table_delimiter
		The character delimiting fields in the output trait table. Valid choices are:tab,space,comma [default: tab]
	`-`-suppress_bifurcating
		If set, don't ensure that tree is fully bifurcating. [default: False]
	-n, `-`-convert_to_nexus
		Convert tree to NEXUS format, including a translate block mapping tip names to numbers. [default: False]
	-c, `-`-convert_values_to_ints
		Convert the values for each character state to integers. [default: False]
	`-`-no_minimum_branch_length
		If set, don't ensure all branches have at least a small but non-zero branchlength. [default: False]
	`-`-supress_tree_filter
		If set, don't filter out tree tips that aren't listed in the trait table [default: False]
	`-`-supress_table_filter
		If set, don't filter out trait table entries that aren't listed in the tree [default: False]
	-r, `-`-add_branch_length_to_root
		Add a short branch to the root node (this is required by some phylogeny programs).  The length of the branch is determined by the min_branch_length option  [default: False]
	-l, `-`-limit_tree_to_otus_fp
		Will prune the reference tree to contain only those tips that are within the given OTU table


**Output:**

Outputs a reformatted tree and trait table.


**Example 1:**

Reformat a tree and trait table with default options:

::

	format_tree_and_trait_table.py -i traits.tab -t tree.nwk -o ./format_output/


