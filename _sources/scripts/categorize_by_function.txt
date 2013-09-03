.. _categorize_by_function:

.. index:: categorize_by_function.py

*categorize_by_function.py* -- Collapse table data to a specified level in a hierarchy.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script collapses hierarchical data to a specified level. For instance, often it is useful to examine KEGG results from a higher level within the pathway hierarchy. Many genes are sometimes involved in multiple pathways, and in these circumstances (also know as a one-to-many relationship), the gene is counted for each pathway. This has a side effect of increasing the total count of genes in the table.


**Usage:** :file:`categorize_by_function.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fp
		The predicted metagenome table
	-o, `-`-output_fp
		The resulting table
	-c, `-`-metadata_category
		The metadata category that describes the hierarchy
	-l, `-`-level
		The level in the hierarchy to collapse to. A value of 0 is not allowed, a value of 1 is the highest level, and any higher value nears the leaves of the hierarchy. For instance, if the hierarchy contains 4 levels, specifying 3 would collapse at one level above being fully specified.
	
	**[OPTIONAL]**
		
	`-`-ignore
		Ignore the comma separated list of names. For instance, specifying --ignore_unknown=unknown,unclassified will ignore those labels while collapsing. The default is to not ignore anything. [default: None]
	-f, `-`-format_tab_delimited
		Output the predicted metagenome table in tab-delimited format [default: False]


**Output:**

Output table is contains gene counts at a higher level within a hierarchy.


Collapse predicted metagenome using KEGG Pathway metadata.

::

	categorize_by_function.py -i predicted_metagenomes.biom -c KEGG_Pathways -l 3 -o predicted_metagenomes.L3.biom

Change output to tab-delimited format (instead of BIOM).

::

	categorize_by_function.py -f -i predicted_metagenomes.biom -c KEGG_Pathways -l 3 -o predicted_metagenomes.L3.txt

Collapse predicted metagenome using taxonomy metadata (not one-to-many).

::

	categorize_by_function.py -i observation_table.biom -c taxonomy -l 1 -o observation_table.L1.biom


