.. _scale_metagenome:

.. index:: scale_metagenome.py

*scale_metagenome.py* -- This script converts metagenomic relative abundance back to sequence counts, by scaling the relative abundnace of each gene in each sample in a biom file by a user-supplied sequencing depth
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`scale_metagenome.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-s, `-`-input_seq_depth_file
		An input tab-delimited table, with samples as the first column and an integer sequencing depth as the second
	-i, `-`-input_count_table
		The input trait counts on per otu basis in biom format (can be gzipped)
	-o, `-`-output_metagenome_table
		The output file for the scaled metagenome


**Output:**

Output is a table of function counts (e.g. KEGG KOs) by sample ids.


Predict metagenomes from genomes.biom and otus.biom.

::

	scale_metagenome.py -i otus.biom -c KEGG_acepic__predict_traits_97.biom.gz -o predicted_metagenomes.biom

Change output format to plain tab-delimited:

::

	scale_metagenome.py -f -i otus.biom -c KEGG_acepic_predict_traits_97.biom.gz -o predicted_metagenomes.tab


