.. _morgan:

======================================
PI-CRUST ala Morgan written in Perl
======================================

This document covers the code in the ``pre_picrust/morgan/`` directory. Questions on usage should be directed to Morgan Langille (morgan.g.i.langille@gmail.com)

Description of code
===================

This code performs the following functions that are the backbone of PICRUST:
 * Creates a reference 16S tree for sequences from reference genomes using PyNAST and RaXML (done once for a given set of genomes; "offline").
 * Performs ancestral state reconstruction (ASR) using various methods and various functional databases for the reference tree (done once per reference genome; "offline").
 * Takes raw metagenomic reads, aligns and filters them using pynast, and places them onto a tree using pplacer. (done once for each metagenome set; "online")
 * Predicts functional abundances for the set of reads in the metagenome sample using a variety of methods including ASR, Nearest Neighbour, and Random.

There is also code to do "leave one out" validation using only reference genomes:
 * Can calculate accuracy (precision and recall) for any of the prediction methods on a genome basis for many functional databases (incl. PFAM, SEED, EC). 
 * Simiarly, calculate accuracy across the various functional classes within a functional database
 * Create numerous scatter plots, violin (distribution) plots, and bar charts. 

Software Requirements
=====================
 * Perl >= v 5.10
 * Perl Modules (BioPerl, Log::Log4perl, JSON, Parallel::ForkManager)
 * PyNAST and GreenGenes template (included with QIIME installation)
 * R
 * R Modules (ape,geiger,vioplot, sciplot)

Installation
============
 * Data resource files are not included with the code, but need to be copied from the shared Dropbox folder. 
"cp genome_test_cases/seed_EC_and_pfam pre_picrust/morgan/data"
 * That's it! You are ready to run PI-CRUST.

Overview
========
 * The main scripts needed to run the pipeline are in the main directory. 
 * These main scripts call other scripts within bin/ which can be accessed but are not well documented. 
 * To get documentation on any of the main scripts you can run the script with the --help option to get the full documentation. (e.g. build_asr.pl --help).
 * Reference trees and other "offline" files will be stored in a sub-directory of the directory "ref_trees". The name of the sub-directory will be named according to the file given to build the reference tree (see usage below for build_ref_16s.pl)
 * Metagenomic read alignments, predictions, and other "online" files will be stored in a sub-directory of "tmp". The name of the sub-directory will be named according to the file containing the metagenomic reads.
 * "run_leave_one_out.pl" is a pipeline script that glues together all of the steps and runs it on a single reference genome leaving that genome out of the reference dataset and treating it's 16S sequence as a single metagenomic read. This script also calls other additional scripts in /bin to test accuracy and create plots. 

Usage examples
==============

This section provides a description of how to run the main code of the PI-CRUST pipeline:

Step 1: Create "offline" files for completed genomes
----------------------------------------------------

Step 1a: Create a reference 16S tree for completed genomes
---------------------------------------------------------

Here's how to create a reference tree from your genome 16S sequences::

    build_ref_16s.pl your_genomes_16S.fa

This takes an input fasta file of 16S sequences from completed genomes. IDs of in the fasta file will be used as genome identifiers. 
Output is written to the directory ``ref_trees/your_genomes_16S.fa``. 
Various ``pynast_*`` and ``RAxML_*`` files are created.


Step 1b: Create ancestral state reconstructions (ASRs) for internal nodes of reference tree
------------------------------------------------------------------------------------------

You'll next need to run ``build_asr.pl`` on the ``your_genomes_16S.fa`` from step 1::

    build_asr.pl -f pfam -m pic -r your_genomes_16S.fa

This performs the users choice of ASR method (pic, ML, REML, etc.) on the reference tree using the user's choice of functional traits (pfam, subsystem, EC, etc.).
Output is written to the directory ``ref_trees/your_genomes_16S.fa``. 
Another version of the reference tree, but with internal node labels added, will created as ``RAxML_result.16s_with_node_labels``.
ASR information will be stored in a compressed tab delimited table (``<func>_<method>_counts.txt.gz``) where column labels correspond to newly created internal tree labels and row labels are function labels. 

Step 2: Make predictions for metagenomic reads
----------------------------------------------

Step 2a: Place metagenomic reads onto reference tree
----------------------------------------------------

The metagenomic reads need to be filtered (by pynast) and placed onto the reference tree (by pplacer)::

    place_reads.pl -r your_genomes_16S.fa meta_reads.fa

This takes the reference tree created in Step 1 and places raw metagenomic 16S reads onto the tree.
Output is written to the directory ``tmp/meta_reads.fa``.
``pynast_trimmed_alignment.jplace`` and ``pplacer.tog.tree`` are json and tree representations,respectively, of the placement output from pplacer.

Step 2b: Make predictions for the metagenomic reads
--------------------------------------------------- 

Make actual functional abundance predictions for each of the functions can be specified::

    make_predictions.pl -f pfam -m pic -r your_genomes_16S.fa -q meta_reads.fa

The command requires the reference tree (created in Step 1 above) and the metagenomic reads name (same as used in step 1a).
The user can choose which functional trait (pfam, subsystem, EC, etc.) and the method (pic, ML, neighbour, random).
Predictions are stored in the directory ``tmp/meta_reads.fa`` in the file ``<func>_<method>_predictions.txt`` (e.g ``pfam_pic_predictions.txt``).

Other Usage
===========

Validate accuracy of PI-CRUST
-----------------------------

Runs the entire validation pipeline including the creation of figures::

    run_leave_one_out.pl -r=your_genomes_16S.fa

Calculates the accuracy of the PI-CRUST software using various functional databases and difference methods for prediction (ASR, NN, random). 
The pipeline removes one genome from the reference genome dataset and treats it as a unknown query (metagenomic) 16S sequence. 
The pipeline makes predictions on this genome and the predicted functional abundances are compared to the known annotated functional abundances for that genome. 
This is repeated for every reference genome in the dataset (>1000 times), thus is very CPU intensive. Some stages such as reference tree construction and ASR may want to be done in pieces on a cluster. Also, if the user has a multi-core machine they can use the --threads option. 

