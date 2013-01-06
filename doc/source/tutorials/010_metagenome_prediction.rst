.. _metagenome_prediction_tutorial:

Metagenome Prediction
=====================

.. include:: ../global.rst


Introduction
------------

This tutorial explains how to predict a microbial community metagenome using PICRUST, based on 16S (or other marker gene) data as detailed in :ref:`otu_picking_tutorial`.
Given this input, the process additionally needs the `PICRUSt software`_ package.

Essential Files
---------------

You should already have the following files before starting this tutorial:

1. The `PICRUSt software`_ package
2. A PICRUSt compatible OTU table (:ref:`otu_picking_tutorial`)

Step 1: Normalize OTU Table
---------------------------

:ref:`normalize_by_copy_number` normalizes the OTU table by dividing each OTU by the known/predicted 16S copy number abdundance.

Input is the users OTU table (that has been referenced picked against Greengenes).

Input and output files are in `biom`_ format::

	normalize_by_copy_number.py 
		-i your_otu_table.biom
		-o your_normalized_otu_table.biom

(Optional) Input format of OTU table can be changed to "classic" `QIIME`_ OTU instead of `biom`_ format using the ``-f`` option: ::

	 normalize_by_copy_number.py 
		-f 
		-i your_otu_table.tab
		-o your_normalized_otu_table.biom

Step 2: Predict Functions For Metagenome
----------------------------------------

:ref:`predict_metagenomes` creates the final metagenome functional predictions. It multiplies each normalized OTU abundance by each predicted functional trait abundance to produce a table of functions (rows) by samples (columns).

Input is the normalized OTU table created by :ref:`normalize_by_copy_number`.

(**NOTE: This step currently requires approx 5GB RAM**)
 
Output is in `biom`_ format by default: ::

	predict_metagenomes.py 
		-i your_normalized_otu_table.biom
		-o your_KEGG_predictions.biom

(Optional) Output format can be changed to tab delimited using ``-f`` option: ::

	predict_metagenomes.py 
		-f 
		-i your_normalized_otu_table.biom
		-o your_KEGG_predictions.tab
