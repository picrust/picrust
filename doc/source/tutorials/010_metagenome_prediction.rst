.. _metagenome_prediction_tutorial:

Metagenome Prediction
=====================

.. include:: ../global.rst


Introduction
------------

This tutorial explains how to predict a microbial community metagenome using PICRUST, based on 16S (or other marker gene) data as detailed in :ref:`otu_picking_tutorial`.
Given this input, the process additionally needs the `PICRUST software`_ package.

Essential Files
---------------

You should already have the following files before starting this tutorial:

1. The `PICRUST software`_ package
2. A PICRUST compatible OTU table (:ref:`otu_picking_tutorial`)
3. The `PICRUST precalculated files`_ **(Note: Not yet available)**

.. warning::

	As of August 2012, our prerelease PICRUST does not yet include precomputed data for all `Greengenes`_ reference OTUs. 
	Therefore you will have to run ``predict_traits.py`` using your OTU table and the ``-l`` option. See **Temporary: Create PICRUST data files** below for more information.

Temporary: Create PICRUST data files
------------------------------------

* Download and unpack `PICRUST temporary files`_.

* If your OTU table is in biom format then you need it also in classic QIIME format.::

	convert_biom.py -b -i your_otu_table.biom -o your_otu_table.tab 

* Now run ``predict_traits`` for 16S copy number: ::

	predict_traits.py
		-t reference_tree.newick
		-l your_otu_table.tab	
 		-i 16S_trait_table.tab	
		-r 16S_asr_wagner_counts.tab
		-o trait_predictions_16S_wagner.biom

* Now run ``predict_traits`` for KEGG: ::

	predict_traits.py 
		-t reference_tree.newick
		-l your_otu_table.tab
		-i KEGG_trait_table.tab
		-r KEGG_asr_wagner_counts.tab
		-o trait_predictions_KEGG_wagner.biom
	


Tutorial Files (optional)
-------------------------

The `PICRUST metagenome tutorial files`_ allow testing of PICRUST using a simple two sample OTU table and smaller support files. Descriptions of the files within this download are as follows:

* OTU table (`biom`_ format ``.biom``)
	* ``hmp_mock_otu_table.biom``
	* This is an OTU table that contains OTU IDs that match IDs in the reference tree.  This is usually done by picking OTUs against the `Greengenes`_ reference package in `QIIME`_.

* Functional trait (`KEGG`_ `KO`_) copy number predictions (`biom`_ format ``.biom``)
    * ``trait_predictions_KEGG_wagner.biom`` is generated using :ref:`predict_traits`
    * Contains predictions for `KO`_ traits for each OTU that has a tip in the reference tree. 

.. warning:

	This file was created specifically for the OTU table included in this tutorial using the ``-l`` option in ``predict_traits.py``.
	It can NOT be used with other OTU tables.

* Marker gene (16S) copy number predictions (`biom`_ format ``.biom``)
    * ``trait_predictions_16S_wagner.biom`` is generated using :ref:`predict_traits`
    * Contains predictions for 16S copy number for each OTU that has a tip in the reference tree. 

.. warning:

	This file was created specifically for the OTU table included in this tutorial using the ``-l`` option in ``predict_traits.py``.
	It can NOT be used with other OTU tables.

Normalize OTU Table
-------------------

:ref:`normalize_by_copy_number` normalizes the OTU table by dividing each OTU by the known/predicted 16S copy number abdundance.

Input is the users OTU table (that has been referenced picked against green genes), and the 16S copy number table for each OTU (produced by :ref:`predict_traits`).

Input and output files are in biom format. ::

	normalize_by_copy_number.py 
		-i your_otu_table.biom
		-c trait_predictions_16S_wagner.biom
		-o your_normalized_otu_table.biom

(Optional) Input format of OTU table can be changed to "classic" `QIIME`_ OTU instead of `biom`_ format using the ``-f`` option: ::

	 normalize_by_copy_number.py 
		-f 
		-i your_otu_table.tab
		-c trait_predictions_16S_wagner.biom
		-o your_normalized_otu_table.biom

Predict Functions For Metagenome
--------------------------------

:ref:`predict_metagenomes` creates the final metagenome functional predictions. It multiplies each normalized OTU abundance by each predicted functional trait abundance to produce a table of functions (rows) by samples (columns).

Input is the normalized OTU table created by :ref:`normalize_by_copy_number` and the functional trait predictions produced by :ref:`predict_traits`. 
 
Output is in `biom`_ format by default: ::

	predict_metagenomes.py 
		-i your_normalized_otu_table.biom
		-c trait_predictions_KEGG_wagner.biom
		-o your_KEGG_predictions.biom

(Optional) Output format can be changed to tab delimited using ``-f`` option: ::

	predict_metagenomes.py 
		-f 
		-i your_normalized_otu_table.biom
		-c trait_predictions_KEGG_wagner.biom
		-o your_KEGG_predictions.tab
