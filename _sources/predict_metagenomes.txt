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
		
	-t, `-`-type_of_prediction
		Type of functional predictions. Valid choices are: ko, cog, rfam [default: ko]
	-g, `-`-gg_version
		Version of GreenGenes that was used for OTU picking. Valid choices are: 13_5, 18may2012 [default: 13_5]
	-c, `-`-input_count_table
		Precalculated function predictions on per otu basis in biom format (can be gzipped). Note: using this option overrides --type_of_prediction and --gg_version. [default: None]
	-a, `-`-accuracy_metrics
		If provided, calculate accuracy metrics for the predicted metagenome.  NOTE: requires that per-genome accuracy metrics were calculated using `predict_traits.py <./predict_traits.html>`_ during genome prediction (e.g. there are "NSTI" values in the genome .biom file metadata)
	`-`-suppress_subset_loading
		Normally, only counts for OTUs present in the sample are loaded.  If this flag is passed, the full biom table is loaded.  This makes no difference for the analysis, but may result in faster load times (at the cost of more memory usage)
	`-`-load_precalc_file_in_biom
		Instead of loading the precalculated file in tab-delimited format (with otu ids as row ids and traits as columns) load the data in biom format (with otu as SampleIds and traits as ObservationIds) [default: False]
	-f, `-`-format_tab_delimited
		Output the predicted metagenome table in tab-delimited format [default: False]


**Output:**

Output is a table of function counts (e.g. KEGG KOs) by sample ids.


Predict KO abundances for a given OTU table picked against the newest version of GreenGenes.

::

	predict_metagenomes.py -i normalized_otus.biom -o predicted_metagenomes.biom

Change output format to plain tab-delimited:

::

	predict_metagenomes.py -f -i normalized_otus.biom -o predicted_metagenomes.txt

Predict COG abundances for a given OTU table.

::

	predict_metagenomes.py -i normalized_otus.biom -t cog -o cog_predicted_metagenomes.biom

Change the version of GG used to pick OTUs

::

	predict_metagenomes.py -i normalized_otus.biom -g 18may2012 -o predicted_metagenomes.biom


