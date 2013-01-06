.. include:: ../global.rst
.. _quickstart_guide:

Quickstart Guide
================

This guide gives a high-level overview of how to use PICRUSt. Before starting, you should have installed PICRUSt as described in :ref:`install`. 

PICRUSt primarily consists of two workflows: gene content inference (detailed in :ref:`genome_prediction_tutorial`) and metagenome inference (detailed in :ref:`metagenome_prediction_tutorial`). Users working with 16S data can use pre-computed gene content information, and as a result don't need to be concerned with the gene content inference workflow. The following sections describe these workflows. Each section links to detailed tutorials that illustrate the exact commands that can be applied as well as example data that you can use to test and learn PICRUSt.

Gene content inference (precomputed for 16S rRNA)
-------------------------------------------------

This workflow infers gene content for OTUs with unknown gene content from OTUs with known content and a phylogenetic tree relating OTUs with known and unknown gene content.

Input:

	* A reference OTU tree (newick format)
	* A gene table for OTUs with known gene composition (i.e., counts of functional genes on a per-OTU basis; biom format by default)

Output:

	* A gene table for OTUs with known and unknown gene composition (i.e., counts of functional genes on a per-OTU basis; biom format by default)

Full details on this workflow, including example data and associated commands for running this workflow can be found in :ref:`genome_prediction_tutorial`.


Metagenome inference
--------------------

This workflow infers gene content for samples from an OTU table and a table of gene contents on a per-OTU basis. The input OTU table must be a closed-reference OTU table, where OTU identifiers correspond to the tips in the reference OTU tree used in the *Gene content inference* workflow. For help with building this OTU table, you should see :ref:`otu_picking_tutorial` which makes use of the `PICRUSt GG reference data`_.

Input: 

	* An OTU table (i.e., counts of OTUs on a per-sample basis; default; biom format by default)

Output:

	* A gene table (i.e., counts of functional genes on a per-sample basis; biom format by default)

Full details on this workflow, including example data and associated commands for running this workflow can be found in :ref:`metagenome_prediction_tutorial`.
