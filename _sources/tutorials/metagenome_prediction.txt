.. _metagenome_prediction_tutorial:

Metagenome Prediction Tutorial
==============================

.. include:: ../global.rst


Introduction
------------

This tutorial explains how to predict a microbial community metagenome using PICRUST, based on 16S (or other marker gene) data as detailed in :ref:`otu_picking_tutorial`.


**BIOM Format**

* Please note that PICRUSt by default uses the relatively new `biom`_ format for representing OTU tables and Gene tables (e.g. KOs by samples). This `has several benefits <http://biom-format.org/documentation/biom_format.html#motivation-for-the-biom-format>`_ including easier integration with other software (e.g. QIIME and others in the future) and allows embedding of extra metadata about both the samples and observations (OTUs/KOs).
* **However, PICRUSt also allow users to input OTU tables and export PICRUSt predictions in tab-delimited format by using the '-f' option (see below).** 
* In addition, users can always convert to/from biom format to tab-delimited format using BIOM's built in `conversion script <http://biom-format.org/documentation/biom_conversion.html>`_.

Requirements
------------

1. You should have already installed PICRUSt (:ref:`install`).
2. A PICRUSt compatible OTU table (:ref:`otu_picking_tutorial`), such as the example file **tutorials/hmp_mock_16S.biom**.

Step 1: Normalize OTU Table
---------------------------

:ref:`normalize_by_copy_number.py <normalize_by_copy_number>` normalizes the OTU table by dividing each OTU by the known/predicted 16S copy number abdundance.

Input is the users OTU table (that has been referenced picked against Greengenes).

Input and output files are in `biom`_ format::

	normalize_by_copy_number.py 
		-i your_otu_table.biom
		-o normalized_otus.biom

**(Optional) Input format of OTU table can be changed to tab-delimited "classic" OTU table instead of BIOM format using the '-f' option:** ::

	 normalize_by_copy_number.py 
		-f 
		-i your_otu_table.tab
		-o normalized_otus.biom

(Optional) Previous examples assume the most recent Greengenes was used for closed OTU picking. Older versions can be specified using the ``--gg_version`` option: ::

          normalize_by_copy_number.py 
	        --gg_version 18may2012
		-i hmp_mock_16S.biom
		-o normalized_otus.biom

Step 2: Predict Functions For Metagenome
----------------------------------------

:ref:`predict_metagenomes.py <predict_metagenomes>` creates the final metagenome functional predictions. It multiplies each normalized OTU abundance by each predicted functional trait abundance to produce a table of functions (rows) by samples (columns).

Input is the normalized OTU table created by :ref:`normalize_by_copy_number.py <normalize_by_copy_number>`.
 
Output is in `biom`_ format by default: ::

	predict_metagenomes.py 
		-i normalized_otus.biom
		-o metagenome_predictions.biom

**(Optional) Output format can be changed from BIOM to tab delimited using '-f' option:** ::

	predict_metagenomes.py 
		-f 
		-i normalized_otus.biom
		-o metagenome_predictions.tab

(Optional) NSTI values for each sample can be obtained using the ``-a`` option: ::

	predict_metagenomes.py 
		-i normalized_otus.biom
		-o metagenome_predictions.tab
		-a nsti_per_sample.tab

(Optional) Previous examples assume the most recent GreenGenes was used for closed OTU picking. Older versions can be specified using the ``--gg_version`` option: ::

        predict_metagenomes.py 
	        --gg_version 18may2012
		-i normalized_otus.biom
		-o metagenome_predictions.biom

(Optional) Previous examples assume that KEGG Orthologs predictions are wanted. Other types of functions (e.g. COGs) can be specified using the ``--type_of_prediction`` option: ::

          predict_metagenomes.py 
	        --type_of_prediction cog
		-i normalized_otus.biom
		-o metagenome_predictions.biom

Step 3: Analysing Predicted Metagenomes
---------------------------------------

Now that you have your predicted metagenome you can analyse it just the same as a real metagenome.
 
Some suggestions for downstream analysis are outlined in :ref:`downstream_analysis_guide`.
