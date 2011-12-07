.. _example:

======================================
Example usage document
======================================

This document covers the code in the ``pre_picrust/caporaso/`` directory. Questions on usage should be directed to Greg Caporaso.

Description of code
===================

This code performs the following functions:
 * Example of how to develop modules
 * Example of how to write unit tests

Usage examples
==============

This section provides a description of how to apply this code to achieve several results.

Step 1: Run the first script...
--------------------------------

Here's how to run the first script::

    script1.py -i input.fasta -o output_dir/

This takes an input fasta file, runs some analyses, and writes the output to ``output_dir``. 


Step 2: Applying the second script to the output of the first...
-----------------------------------------------------------------

You'll next need to run ``script2.py`` on the ``output_dir`` from step 1::

    script2.py -i output_dir/ -o step2_out.txt

This performs various statistical analyses on the output of step 1 and writes the results to ``step2_out.txt``.
 