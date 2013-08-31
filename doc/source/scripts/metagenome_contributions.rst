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
	-o, `-`-output_fp
		The output file for the metagenome contributions
	
	**[OPTIONAL]**
		
	-t, `-`-type_of_prediction
		Type of functional predictions. Valid choices are: ko, cog, rfam [default: ko]
	-g, `-`-gg_version
		Version of GreenGenes that was used for OTU picking. Valid choices are: 13_5, 18may2012 [default: 13_5]
	-c, `-`-input_count_table
		Precalculated function predictions on per otu basis in biom format (can be gzipped). Note: using this option overrides --type_of_prediction and --gg_version. [default: None]
	`-`-suppress_subset_loading
		Normally, only counts for OTUs present in the sample are loaded.  If this flag is passed, the full biom table is loaded.  This makes no difference for the analysis, but may result in faster load times (at the cost of more memory usage)
	`-`-load_precalc_file_in_biom
		Instead of loading the precalculated file in tab-delimited format (with otu ids as row ids and traits as columns) load the data in biom format (with otu as SampleIds and traits as ObservationIds) [default: False]
	-l, `-`-limit_to_function
		If provided, only output predictions for the specified function ids.  Multiple function ids can be passed using comma delimiters.


**Output:**

Output is a tab-delimited column indicating OTU contribution to each function.


Partition the predicted contribution to the  metagenomes from each organism in the given OTU table, limited to only K00001, K00002, and K00004.

::

	metagenome_contributions.py -i normalized_otus.biom -l K00001,K00002,K00004 -o ko_metagenome_contributions.tab

Partition the predicted contribution to the  metagenomes from each organism in the given OTU table, limited to only COG0001 and COG0002.

::

	metagenome_contributions.py -i normalized_otus.biom -l COG0001,COG0002 -t cog -o cog_metagenome_contributions.tab


