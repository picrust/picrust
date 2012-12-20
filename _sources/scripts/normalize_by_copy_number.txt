.. _normalize_by_copy_number:

===========================
normalize_by_copy_number.py
===========================

Description of code
===================

 * The ``normalize_by_copy_number.py`` script normalizes an OTU table by the predicted 16S copy number abundances for each OTU.

Input Requirements
==================
 * An OTU table (that has been properly reference picked against green genes) in biom format. This table should have SampleIds (columns) as some arbitrary sample id, and ObservationIds (rows) as green gene identifiers.
 * A table of pre-computed 16S copy number abundances for all green gene identifiers (in the gg ref package). Usually downloaded from PICRUST (**16S_acepic_predict_traits_97.biom.gz**).Also in biom format.

Output
======
 * A normalized OTU table (biom format).

Mandatory Options
=================
 * ``-i``: the input OTU table
 * ``-c``: the input pre-computed 16S copy number abundances
 * ``-o``: the normalized OTU table

Optional Options
================
 * ``-f``: input OTU table in tab-delimited format (instead of biom).


Usage examples
==============

This section provides a description of how to run ``normalize_by_copy_number.py``:

* Basic usage::

    predict_metagenomes.py -i your_otu_table.biom -c 16S_acepic_predict_traits_97.biom.gz -o your_normalized_otu_table.biom 

* Input format for OTU table can be set to plain tab-delimited format (instead of biom) using ``-f`` option::

    predict_metagenomes.py -f -i your_otu_table.tab -c KEGG_acepic_predict_traits_97.biom.gz -o your_normalized_otu_table.biom 

* ``-c`` file does can be gzipped or uncompressed::

    predict_metagenomes.py -i your_otu_table.biom -c KEGG_acepic_predict_traits_97.biom -o your_normalized_otu_table.biom 
