.. _metagenome_contributions:

.. index:: metagenome_contributions.py

*metagenome_contributions.py* -- This script partitions metagenome functional contributions according to function, OTU, and sample, for a given OTU table.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`metagenome_contributions.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_otu_table
		The input otu table in biom format
	-c, `-`-input_count_table
		The input trait counts on per otu basis in biom format (can be gzipped)
	-o, `-`-output_metagenome_table
		The output file for the predicted metagenome
	
	**[OPTIONAL]**
		
	-a, `-`-accuracy_metrics
		If provided, calculate accuracy metrics for the predicted metagenome.  NOTE: requires that per-genome accuracy metrics were calculated using `predict_traits.py <./predict_traits.html>`_ during genome prediction (e.g. there are "NSTI" values in the genome .biom file metadata)
	`-`-limit_to_function
		If provided, only output predictions for the specified function ids.  Multiple function ids can be passed using comma delimiters.
	-f, `-`-format_tab_delimited
		Output the predicted metagenome table in tab-delimited format [default: False]


**Output:**

Output is a table of function counts (e.g. KEGG KOs) by sample ids.


Partition the predicted contribution to the  metagenomes from each organism in otus.biom, using the predicted genes for each organism in genes.biom.

::

	metagenome_contributions.py -i otus.biom -c KEGG_acepic__predict_traits_97.biom.gz -o predicted_metagenomes.biom

Change output format to plain tab-delimited:

::

	metagenome_contributions.py -f -i otus.biom -c KEGG_acepic_predict_traits_97.biom.gz -o predicted_metagenomes.tab


