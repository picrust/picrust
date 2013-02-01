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


**Output:**

Output table is contains gene counts at a higher level within a hierarchy.


Collapse predicted metagenome results.

::

	categorize_by_function.py -i metagenome.biom -c "KEGG Pathways" -l 3 -o metagenome_at_level3.biom


