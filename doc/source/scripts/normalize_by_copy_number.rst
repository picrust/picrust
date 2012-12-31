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
 * A table of pre-computed 16S copy number abundances for all green gene identifiers (in the gg ref package).  This file is generated using ``predict_traits.py``, but for most use-cases is precalculated ('data/16S_precalculated.biom.gz') and used by default.

Output
======
 * A normalized OTU table (biom format).

Mandatory Options
=================
 * ``-i``: the input OTU table
 * ``-o``: the normalized OTU table

Optional Options
================
 * ``-f``: input OTU table in tab-delimited format (instead of biom).
 * ``-c``: the input pre-computed 16S copy number abundances


Usage examples
==============

This section provides a description of how to run ``normalize_by_copy_number.py``:

* Basic usage::

    normalize_by_copy_number.py -i your_otu_table.biom -o normalized_otus.biom 

* Input format for OTU table can be set to plain tab-delimited format (instead of biom) using ``-f`` option::

    normalize_by_copy_number.py -f -i your_otu_table.tab -o normalized_otus.biom 

* A different set of marker gene copy numbers can be specified using ``-c`` and this file can be gzipped or uncompressed::

    normalize_by_copy_number.py -i your_otu_table.biom -c your_marker_gene_copy_number.biom -o normalized_otus.biom 
