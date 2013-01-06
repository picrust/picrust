.. _predict_metagenomes:

.. index:: predict_metagenomes.py

*predict_metagenomes.py* -- This script produces the actual metagenome functional predictions for a given OTU table.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`predict_metagenomes.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_otu_table
		The input otu table in biom format
	-o, `-`-output_metagenome_table
		The output file for the predicted metagenome
	
	**[OPTIONAL]**
		
	-c, `-`-input_count_table
		Precalculated function predictions on per otu basis in biom format (can be gzipped) [default: /Users/caporaso/.virtualenvs/picrust/lib/python2.7/site-packages/picrust/data/ko_precalculated.biom.gz]
	-a, `-`-accuracy_metrics
		If provided, calculate accuracy metrics for the predicted metagenome.  NOTE: requires that per-genome accuracy metrics were calculated using `predict_traits.py <./predict_traits.html>`_ during genome prediction (e.g. there are "NSTI" values in the genome .biom file metadata)
	-f, `-`-format_tab_delimited
		Output the predicted metagenome table in tab-delimited format [default: False]


**Output:**

Output is a table of function counts (e.g. KEGG KOs) by sample ids.


Predict metagenomes from genomes.biom and otus.biom.

::

	predict_metagenomes.py -i normalized_otus.biom -o predicted_metagenomes.biom

Change output format to plain tab-delimited:

::

	predict_metagenomes.py -f -i normalized_otus.biom -o predicted_metagenomes.tab


