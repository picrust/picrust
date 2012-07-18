.. _metagenome_prediction_tutorial:

PICRUST Metagenome Prediction Tutorial
======================================

Introduction
------------
This tutorial explains the final steps in the **PI-CRUST** pipeline that takes in a users OTU table and produces metagenome predictions. 

These commands would be done after the steps outlined in the `PICRUST Genome Prediction Tutorial <./genome_prediction.html>`_. 

Essential Files
---------------

In this tutorial, you will need an OTU table and the output files from the 'predict_traits.py' script. 

The steps required to get the output from the intial steps in the PI-CRUST pipeline are outlined in the :ref:`genome_prediction_tutorial`.

**Note: If you are trying to run PICRUST with your own OTU table ensure you have properly reference picked your OTU table:** :ref:`otu_picking_tutorial`

**Note2: Trait predictions have not been pre-calculated yet for all greengenes tips. Therefore you will have to run 'predict_traits.py' using your OTU table using the `-l` option. See:** :ref:`genome_prediction_tutorial`

You can download all files required to run this tutorial here (https://www.dropbox.com/s/yndl609h5m069y6/metagenome_prediction_tutorial_files.zip). Descriptions of these files are below. 

* Marker gene (16S) copy number predictions (biom format .biom)
    * Contains predictions for 16S copy number for each OTU that has a tip in the reference tree. 
    * We will use the precomputed file **trait_predictions_16S_wagner.biom**
    * *Note: This file was created specifically for the OTU table included in this tutorial using the '-l' option in 'predict_traits.py'. It can NOT be used with other OTU tables.* 

* Functional trait (KEGG KO) copy number predictions (biom format .biom)
    * Contains predictions for KO traits for each OTU that has a tip in the reference tree. 
    * We will use the precomputed file **trait_predictions_KEGG_wagner.biom**.
    * *Note: This file was created specifically for the OTU table included in this tutorial using the '-l' option in 'predict_traits.py'. It can NOT be used with other OTU tables.*

* OTU table (biom format .biom)
    * This is an OTU table that contains OTU ids that match ids in the reference tree. This is usally done by picking OTUs against the GG reference package in QIIME.
    * Example OTU table  **hmp_mock_otu_table.biom**.


Normalize OTU Table
-------------------
`normalize_by_copy_number.py <../scripts/normalize_by_copy_number.html>`_ normalizes the OTU table by dividing each OTU by the known/predicted 16S copy number abdundance.

Input is the users OTU table (that has been referenced picked against green genes), and the 16S copy number table for each OTU (produced by 'predict_traits.py').

Input and output files are in biom format. ::

	normalize_by_copy_number.py -i hmp_mock_otu_table.biom -c trait_predictions_16S_wagner.biom -o normalized/hmp_mock_wagner.biom

(Optional) Input format of OTU table can be changed to "classic" QIIME OTU instead of biom format using the '-f' option: ::

	 normalize_by_copy_number.py -f -i hmp_mock_otu_table.biom -c trait_predictions_16S_wagner.biom -o normalized/hmp_mock_wagner.biom

Predict Functions For Metagenome
--------------------------------
`predict_metagenomes.py <../scripts/predict_metagenomes.html>`_ creates the final metagenome functional predictions. It multiplies each normalized OTU abundance by each predicted functional trait abundance to produce a table of functions (rows) by samples (columns).

Input is the normalized OTU table created by 'normalize_by_copy_number.py' and the functional trait predictions produced by 'predict_traits.py'. 
 
Output is in 'biom' format by default: ::

	predict_metagenomes.py -i normalized/hmp_mock_wagner.biom -g trait_predictions_KEGG_wagner.biom -o hmp_mock_predictions_wagner.biom

(Optional) Output format can be changed to tab delimited using '-f' option: ::

	predict_metagenomes.py -f -i normalized/hmp_mock_wagner.biom -g trait_predictions_KEGG_wagner.biom -o hmp_mock_predictions_wagner.tab
