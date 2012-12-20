.. _dan_and_jesse:

===============================================================
R-based scripts for genome reconstruction using R package "ace"
===============================================================

This document covers the code in the ``pre_picrust/dan_and_jesse/r_scripts`` directory. Questions on usage should be directed to Dan Knights (danknights@gmail.com)

Description of code
===================

The file "predict_tip_states.r" implements a wrapper for phylogeny-informed organism trait reconstruction via ancestral state reconstruction, k-nearest neighbors. All predictions are made by removing one organism from the tree and predicting its states using all remaining organisms. Ancestral state reconstruction is performed using the 'ace' function from the R library "ape".

Software Requirements
=====================
 * R
 * "ape" package in R (install by running the command "install.packages('ape')" inside R)
 * "Matrix" package in R

Usage examples
==============

Input files:
 * Newick-based tree file, must contain branch lengths
 * Tab-delimited organism state table. Rows are organims, columns are traits, 0 indicates "absent", 1 indicates "present".

Output files:
 * Tab-delimited table of leave-one-out predictions for organism trait states as probabilities
 * Inferred transition rate matrix when using "ape".

Get help::

	R --slave --vanilla --args -h < predict_tip_states.r 

Leave-one-out predictions using ancestral state reconstruction::

	R --slave --vanilla --args -t test_data/randomtree.newick -s test_data/randomstates.txt -m ace -o res < predict_tip_states.r 

Leave-one-out predictions using with ancestral state reconstruction with fixed transition rate::
	R --slave --vanilla --args -t test_data/randomtree.newick -s test_data/randomstates.txt -m ace -o res --fixed_rate .1 < predict_tip_states.r 

Leave-one-out predictions using k-nearest neighbors with k=1::

	R --slave --vanilla --args -t test_data/randomtree.newick -s test_data/randomstates.txt -m knn -o res < predict_tip_states.r 

Leave-one-out predictions using k-nearest neighbors with k=3::

	R --slave --vanilla --args -t test_data/randomtree.newick -s test_data/randomstates.txt -m knn -k 3 -o res < predict_tip_states.r 
