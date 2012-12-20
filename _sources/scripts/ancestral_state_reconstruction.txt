.. _ancestral_state_reconstruction:

==============================
Ancestral State Reconstruction
==============================

This document covers use of the ``picrust/ancestral_state_reconstruction.py`` script with some information about the various app controllers that it calls. Questions on usage should be directed to Morgan Langille (morgan.g.i.langille@gmail.com).

Description of code
===================

 * The ``ancestral_state_reconstruction.py`` script is a wrapper script that interfaces with various third party programs of ancestral state recontruction (ASR).
 * Each of the third party ASR programs all have their own pycogent app controllers. 
 * Each of these app controllers have a method called ``XXXXX_for_picrust``, where ``XXXXX`` is the name of ASR method. This is called by ``ancestral_state_reconstruction.py`` and ensures that same input and output formats as agreed to by the PICRUST developers.

Input Requirements
==================
 * A tree in newick format with unique labels for the tips AND the internal nodes. The tree must be dichotomous (e.g. no polytomies), with no zero branch lengths. 
 * A tab-delimited table with the first row being a list of trait identifiers, the first column being a list of ALL tip labels from the tree, and the cell values being counts of data.
 * Note: These files will usually be the outputs produced by the ``format_tree_and_trait_table.py`` script. 

Output
======
 * A tab-delmited table with the first row being a list of trait identifiers, the second column being a list of all internal node labels from the tree, and the cell values being the predicted counts of data. 
 * If a the ASR method provides confidence information (posterior probabilities, confidence intervals, etc.), then an additional table will be output that provides the same columns and row headings, but the cell values will contain probability information for the predictions. 

Mandatory Options
=================
 * ``-i``: the input trait table described above
 * ``-t``: the input tree described above

Optional Options
================
 * ``-m``: the asr method to be used. Can be one of (``ace_ml``,``ace_reml``, ``ace_pic``, ``wagner``). Default is ``wagner``.
 * ``-o``: name for output table containing ASR predictions of just the counts [default:asr_counts.tab]


Usage examples
==============

This section provides a description of how to run ``ancestral_state_reconstruction.py``:

Basic Example::

    ancestral_state_reconstruction.py -i your_trait_table.txt -t your_tree.newick

Specify a different ASR method::

    ancestral_state_reconstruction.py -i your_trait_table.txt -t your_tree.newick -m ace_pic

Specify a different a name for the output file(s)::

    ancestral_state_reconstruction.py -i your_trait_table.txt -t your_tree.newick -o your_output_name
    
