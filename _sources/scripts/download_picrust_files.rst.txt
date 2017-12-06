.. _download_picrust_files:

.. index:: download_picrust_files.py

*download_picrust_files.py* -- Download PICRUSt precalculated files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Downloads the PICRUSt precalculated files and stores them in the PICRUSt data directory.

Note that the script will skip downloading files already in the data directory (unless the ``--force`` argument is used).

**Usage:** :file:`download_picrust_files.py [options]`

**Input Arguments:**

.. note::

	**[OPTIONAL]**

	`-`-force
		Download the files again (i.e. overwrite any existing ones) [default: False]
	-g, `-`-gg_version
		Version of GreenGenes that was used for OTU picking. Valid choices are: 13_5, 18may2012 [default: 13_5]
	-t, `-`-type_of_prediction
		Type of functional predictions. Valid choices are: ko, cog, rfam [default: ko]
	`-`-with_confidence
		Download the precalculated files that allow confidence interval calculations for metagenome predictions. Only available for GreenGenes 13_5. [default: False]

**Output:**

Files are downloaded and stored in PICRUSt's internally used data directory.

Download the default files (GreenGenes 13_5, KEGG KOs):

::

	download_picrust_files.py

Download the default files (GreenGenes 13_5, KEGG KOs) *and* the files necessary for calculating confidence intervals:

::

	download_picrust_files.py --with_confidence

Use GreenGenes 18may2012 and COG:

::

	download_picrust_files.py -g 18may2012 -t cog
