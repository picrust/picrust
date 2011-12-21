.. _josh:

=====================================
Pre-PI-CRUST from the Huttenhower Lab
=====================================

This document covers the code in the ``pre_picrust/huttenhower/`` directory. Questions on usage should be directed to Joshua Reyes (jreyes@post.harvard.edu).

Description of code
===================

This code performs the following functions that are the backbone of PI-CRUST:
 * Constructs a (NCBI) labeled phylogeny from the input OTU abundance table.
 * Performs ancestral state reconstruction from the tree using a decaying exponential weighted vote to produce a gene count for each missing genome.
 * Predicts functional abundances simply by multiplying OTU abundance by gene count and summing across all OTUs.

Software Requirements
=====================
 * Python >= 2.6
 * Numpy, a package for scientific computing in Python
 * ETE 2, a package for the manipulation, analysis, and visualization of phylogenetic trees
 * SCons, a software construction tool similar to the Make utility
 
Installation
============
 To install, make ``$PICRUST_HOME/bin/ANFunc`` path accessible. 
 
 That's it! You are ready to run PI-CRUST.

Overview
========
 * To run the pipeline simply copy as many OTU abundance PCL files as you like into the ``input`` directory and let SCons take over::

    scons

 * Final results are saved in ``output/<input name>_kos.pcl`` by default.

Notes
=====
 * The SConstruct file drives the entire pipeline and can be tweaked to suit your needs. 
 * This pipeline relies on a list of taxa (``input/taxdump.txt``) and a count of KOs for each genome from NCBI (``genomes_02-zro-flt.txt``) and are included. 
 * Intermediate calculations are saved in ``output`` folder. Trees are saved using an internal representation in ``.fir`` files. Intermediate tables have end in the suffix ``_fir.pcl``.
