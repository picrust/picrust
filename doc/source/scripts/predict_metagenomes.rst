.. _predict_metagenomes:

======================
predict_metagenomes.py
======================

This document covers use of the ``predict_metagenomes.py`` script.

Description of code
===================

 * The ``predict_metagenomes.py`` script produces the actual metagenome functional predictions for a given OTU table.

Input Requirements
==================
 * An OTU table (that has been properly reference picked against green genes) in biom format. This table should have SampleIds (columns) as some arbitrary sample id, and ObservationIds (rows) as green gene identifiers.
 * A table of pre-computed genome contents for all green gene identifiers (in the gg ref package). Usually downloaded from PICRUST (**KEGG_acepic_predict_traits_97.biom.gz**).Also in biom format.

Output
======
 * A table of functional identifiers (e.g. KEGG KOs) by the users sample ids. (biom format).

Mandatory Options
=================
 * ``-i``: the input OTU table
 * ``-c``: the input pre-computed genome contents
 * ``-o``: the output file

Optional Options
================
 * ``-f``: output the predicted metagenome table in tab-delimited format (instead of biom).


Usage examples
==============

This section provides a description of how to run ``predict_metagenomes.py``:

* Basic usage::

    predict_metagenomes.py -i your_normalized_otu_table.biom -c KEGG_acepic_predict_traits_97.biom.gz -o your_output.biom 

* Output format can be set to plain tab-delimited format (instead of biom) using ``-f`` option::

    predict_metagenomes.py -f -i your_normalized_otu_table.biom -c KEGG_acepic_predict_traits_97.biom.gz -o your_output.tab 

* ``-c`` file does can be gzipped or uncompressed::

    predict_metagenomes.py -i your_normalized_otu_table.biom -c KEGG_acepic_predict_traits_97.biom -o your_output.biom 
    
