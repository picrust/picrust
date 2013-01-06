.. _normalize_by_copy_number:

.. index:: normalize_by_copy_number.py

*normalize_by_copy_number.py* -- Normalize an OTU table by marker gene copy number
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`normalize_by_copy_number.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_otu_fp
		The input otu table filepath in biom format
	-o, `-`-output_otu_fp
		The output otu table filepath in biom format
	
	**[OPTIONAL]**
		
	-c, `-`-input_count_fp
		The input marker gene counts on per otu basis in biom format (can be gzipped) [default: /Users/caporaso/.virtualenvs/picrust/lib/python2.7/site-packages/picrust/data/16S_precalculated.biom.gz]
	`-`-metadata_identifer
		Identifier for copy number entry as observation metadata [default: CopyNumber]
	-f, `-`-input_format_classic
		Input otu table (--input_otu_fp) is in classic Qiime format [default: False]


**Output:**

A normalized OTU table


Normalize the counts in raw_otus.biom. Write the resulting table to normalized_otus.biom.

::

	normalize_by_copy_number.py -i raw_otus.biom -o normalized_otus.biom

Input tab-delimited OTU table:

::

	normalize_by_copy_number.py -f -i raw_otus.tab -o predicted_metagenomes.biom


