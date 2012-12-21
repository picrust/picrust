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
 * A table of pre-computed genome contents for all green gene identifiers (in the gg ref package). This file is generated using ``predict_traits.py``, but for most use-cases is precalculated ('data/ko_precalculated.biom.gz') and used by default.

Output
======
 * A table of functional identifiers (e.g. KEGG KOs) by the users sample ids. (biom format).

Mandatory Options
=================
 * ``-i``: the input OTU table
 * ``-o``: the output file

Optional Options
================
 * ``-f``: output the predicted metagenome table in tab-delimited format (instead of biom).
 * ``-c``: the input pre-computed genome contents (biom format and optionally gzipped).

Usage examples
==============

This section provides a description of how to run ``predict_metagenomes.py``:

* Basic usage::

    predict_metagenomes.py -i normalized_otus.biom -o predicted_metagenomes.biom 

* Output format can be set to plain tab-delimited format (instead of biom) using ``-f`` option::

    predict_metagenomes.py -f -i normalized_otus.biom -o predicted_metagenomes.biom 

* A different set of functions can be specified using ``-c`` and this file can be gzipped or uncompressed::

    predict_metagenomes.py -i normalized_otus.biom -c your_output_from_predict_traits.biom -o predicted_metagenomes.biom 
    
