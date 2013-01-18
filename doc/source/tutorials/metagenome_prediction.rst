.. _metagenome_prediction_tutorial:

Metagenome Prediction Tutorial
==============================

.. include:: ../global.rst


Introduction
------------

This tutorial explains how to predict a microbial community metagenome using PICRUST, based on 16S (or other marker gene) data as detailed in :ref:`otu_picking_tutorial`.

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
		-i hmp_mock_16S.biom
		-o normalized_otus.biom

(Optional) Input format of OTU table can be changed to "classic" `QIIME`_ OTU instead of `biom`_ format using the ``-f`` option: ::

	 normalize_by_copy_number.py 
		-f 
		-i hmp_mock_16S.tab
		-o normalized_otus.biom

Step 2: Predict Functions For Metagenome
----------------------------------------

:ref:`predict_metagenomes.py <predict_metagenomes>` creates the final metagenome functional predictions. It multiplies each normalized OTU abundance by each predicted functional trait abundance to produce a table of functions (rows) by samples (columns).

Input is the normalized OTU table created by :ref:`normalize_by_copy_number.py <normalize_by_copy_number>`.

(**NOTE: This step currently requires approx 5GB RAM**)
 
Output is in `biom`_ format by default: ::

	predict_metagenomes.py 
		-i normalized_otus.biom
		-o metagenome_predictions.biom

(Optional) Output format can be changed to tab delimited using ``-f`` option: ::

	predict_metagenomes.py 
		-f 
		-i normalized_otus.biom
		-o metagenome_predictions.tab

(Optional) NSTI values for each sample can be obtained using the ``-a`` option: ::

	predict_metagenomes.py 
		-i normalized_otus.biom
		-o metagenome_predictions.tab
		-a nsti_per_sample.tab

Step 3: Analysing Predicted Metagenomes
---------------------------------------

Know that you have your predicted metagenome you can analyse it just the same as a real metagenome. 
Some suggestions for downstream analysis are outlined in :ref:`downstream_analysis_guide`.
