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
		
	-g, `-`-gg_version
		Version of GreenGenes that was used for OTU picking. Valid choices are: 18may2012 [default: 18may2012]
	-c, `-`-input_count_fp
		Precalculated input marker gene copy number predictions on per otu basis in biom format (can be gzipped).Note: using this option overrides --gg_version. [default: None]
	`-`-metadata_identifer
		Identifier for copy number entry as observation metadata [default: CopyNumber]
	-f, `-`-input_format_classic
		Input otu table (--input_otu_fp) is in classic Qiime format [default: False]
	`-`-load_precalc_file_in_biom
		Instead of loading the precalculated file in tab-delimited format (with otu ids as row ids and traits as columns) load the data in biom format (with otu as SampleIds and traits as ObservationIds) [default: False]


**Output:**

A normalized OTU table


Normalize the counts in raw_otus.biom. Write the resulting table to normalized_otus_from_biom.biom.

::

	normalize_by_copy_number.py -i hmp_mock_16S.biom -o normalized_otus_from_biom.biom

Input tab-delimited OTU table:

::

	normalize_by_copy_number.py -f -i hmp_mock_16S.txt -o normalized_otus_from_txt.biom


