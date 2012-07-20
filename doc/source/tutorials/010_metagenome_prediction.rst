.. _metagenome_prediction_tutorial:

.. include:: ../global.rst

Metagenome Prediction
=====================

Introduction
------------

This tutorial explains how to predict a microbial community metagenome using PICRUST, based on 16S (or other marker gene) data as detailed in :ref:`otu_picking_tutorial`.
Given this input, the process additionally needs the `PICRUST software`_ package.

Essential Files
---------------

In this tutorial, you will need an OTU table and the output files from the ``predict_traits.py`` script. 

The steps required to get the output from the intial steps in the PICRUST pipeline are outlined in :ref:`genome_prediction_tutorial`.

.. warning::

	If you are trying to run PICRUST with your own OTU table, ensure you have properly reference picked your OTU table: :ref:`otu_picking_tutorial`

.. warning::

	As of July 2012, our prerelease PICRUST does not yet include precomputed data for all `Greengenes`_ reference OTUs. 
	Therefore you will have to run ``predict_traits.py`` using your OTU table and the ``-l`` option. See: :ref:`genome_prediction_tutorial`

You can download all files required to run this tutorial here: https://www.dropbox.com/s/yndl609h5m069y6/metagenome_prediction_tutorial_files.zip.
Descriptions of these files are below. 

* OTU table (`biom`_ format ``.biom``)
	* ``hmp_mock_otu_table.biom``
	* This is an OTU table that contains OTU IDs that match IDs in the reference tree.  This is usually done by picking OTUs against the `Greengenes`_ reference package in `QIIME`_.

* Functional trait (`KEGG`_ `KO`_) copy number predictions (`biom`_ format ``.biom``)
    * ``trait_predictions_KEGG_wagner.biom`` is precomputed by PICRUST.
    * Contains predictions for `KO`_ traits for each OTU that has a tip in the reference tree. 

.. warning:

	This file was created specifically for the OTU table included in this tutorial using the ``-l`` option in ``predict_traits.py``.
	It can NOT be used with other OTU tables.

* Marker gene (16S) copy number predictions (`biom`_ format ``.biom``)
    * ``trait_predictions_16S_wagner.biom``
    * Contains predictions for 16S copy number for each OTU that has a tip in the reference tree. 

.. warning:

	This file was created specifically for the OTU table included in this tutorial using the ``-l`` option in ``predict_traits.py``.
	It can NOT be used with other OTU tables.

Normalize OTU Table
-------------------

:ref:`normalize_by_copy_number` normalizes the OTU table by dividing each OTU by the known/predicted 16S copy number abdundance.

Input is the users OTU table (that has been referenced picked against green genes), and the 16S copy number table for each OTU (produced by :ref:`predict_traits`).

Input and output files are in biom format. ::

	normalize_by_copy_number.py -i hmp_mock_otu_table.biom
		-c trait_predictions_16S_wagner.biom
		-o normalized/hmp_mock_wagner.biom

(Optional) Input format of OTU table can be changed to "classic" `QIIME`_ OTU instead of `biom`_ format using the ``-f`` option: ::

	 normalize_by_copy_number.py -f -i hmp_mock_otu_table.biom
		-c trait_predictions_16S_wagner.biom
		-o normalized/hmp_mock_wagner.biom

Predict Functions For Metagenome
--------------------------------

:ref:`predict_metagenomes` creates the final metagenome functional predictions. It multiplies each normalized OTU abundance by each predicted functional trait abundance to produce a table of functions (rows) by samples (columns).

Input is the normalized OTU table created by :ref:`normalize_by_copy_number` and the functional trait predictions produced by :ref:`predict_traits`. 
 
Output is in `biom`_ format by default: ::

	predict_metagenomes.py -i normalized/hmp_mock_wagner.biom
		-g trait_predictions_KEGG_wagner.biom
		-o hmp_mock_predictions_wagner.biom

(Optional) Output format can be changed to tab delimited using ``-f`` option: ::

	predict_metagenomes.py -f -i normalized/hmp_mock_wagner.biom
		-g trait_predictions_KEGG_wagner.biom
		-o hmp_mock_predictions_wagner.tab
