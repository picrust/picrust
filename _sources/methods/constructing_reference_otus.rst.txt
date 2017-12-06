.. _constructing_reference_otus:

Constructing PICRUSt reference OTU collection
==============================================

Introduction
------------
This document covers how the nested reference OTU collection was constructed for PICRUST. The script and associated test code are available in the ``qiime-dev/nested_reference_otus`` git repository located `here <https://github.com/qiime-dev/nested_reference_otus>`_. Download links for these reference collections are:

 * img_gg_otus_18may2012.tgz (`download <https://s3.amazonaws.com/picrust-public-data/img_gg_otus_18may2012.tgz>`_ ; md5: 5f1de69a310dc0c261ea76b48ba5dbe2)

Overall workflow design
-----------------------
Beginning with a collection of sequences of a single marker gene (16S in our case), a user-defined list of percent identities to cluster at, and an optional phylogenetic tree, the following steps are performed:

 #. Pick OTUs (using ``pick_otus.py`` in `QIIME <http://www.qiime.org>`_) on input sequences at first specified percent identity using `uclust <http://www.drive5.com/uclust>`_, suppressing QIIME's pre-sorting of sequences (by passing ``-DB`` to ``pick_otus.py``). Suppressing of the pre-sorting is included so the user can pre-sort sequences as they prefer before OTU picking: sequences that are listed first in the users input will preferentially be cluster centroids. In the case of PICRUST, these are sequences for which we have full genomes.
 #. Pick a representative sequence for each OTU as the centroid of the OTU cluster. 
 #. Filter the input tree (if any) to contain only tips that are representative sequences generated in the previous step. If no input tree was provided, do not create a tree.
 #. Using the representative sequence collection from the previous step, pick otus on input sequences at the next specified percent identity using the same parameters as applied in Step1. 
 #. Pick a representative sequence for each OTU as the centroid of the OTU cluster.
 #. Filter the tree generated in the previous iteration (if any) to contain only tips that correspond to representative sequences in the previous step. If no input tree was generated in the previous step, do not create a tree.
 #. Return to step 4 until OTUs have been picked at all specified percent identities.

The output of this workflow will be OTU representative sequences and corresponding phylogenetic trees (if applicable) at the user-specified percent identities. The OTUs will be nested such that sequences which are clustered in the same OTU at a given percent identity are guaranteed to be clustered in the same OTU at all higher percent identities.

Notes on specific versions of reference OTU builds
--------------------------------------------------

img_gg_otus_18may2012
^^^^^^^^^^^^^^^^^^^^^

 * Run with QIIME 1.5.0 on ami-e4bf1b8d
 * Run by Greg Caporaso on 18 May 2012
 * nested_reference_workflow.py commit 73799024fa4dc7d26d33743ddcb4c1c624520e23
 * Input MD5 sums:
    * IMG_v350_sel4cni.lanemask.min1250nt.rooted.t2tGGrelease.ntree : e988ef28ce22cfa2adab742b1bf0ad98
    * IMG_v350_sel4cni_unaligned.fasta : d3ac90d3b6f24685cac056d1907daa7c
 * The full tree contained single quotes around the IMG sequence identifiers, so I had to commit a fix to QIIME to strip these from the tip names before filtering. I then re-ran filtering using QIIME 1.5.0-dev revision 3062 on my laptop::

	filter_tree.py -i IMG_v350_sel4cni.lanemask.min1250nt.rooted.t2tGGrelease.ntree -f img_gg_otus_18may2012/rep_set/99_otus_img_gg_18may2012.fasta -o img_gg_otus_18may2012/trees/99_otus_img_gg_18may2012.tre
	filter_tree.py -i img_gg_otus_18may2012/trees/99_otus_img_gg_18may2012.tre -f img_gg_otus_18may2012/rep_set/97_otus_img_gg_18may2012.fasta -o img_gg_otus_18may2012/trees/97_otus_img_gg_18may2012.tre
	filter_tree.py -i img_gg_otus_18may2012/trees/97_otus_img_gg_18may2012.tre -f img_gg_otus_18may2012/rep_set/94_otus_img_gg_18may2012.fasta -o img_gg_otus_18may2012/trees/94_otus_img_gg_18may2012.tre
	filter_tree.py -i img_gg_otus_18may2012/trees/94_otus_img_gg_18may2012.tre -f img_gg_otus_18may2012/rep_set/91_otus_img_gg_18may2012.fasta -o img_gg_otus_18may2012/trees/91_otus_img_gg_18may2012.tre
	filter_tree.py -i img_gg_otus_18may2012/trees/91_otus_img_gg_18may2012.tre -f img_gg_otus_18may2012/rep_set/88_otus_img_gg_18may2012.fasta -o img_gg_otus_18may2012/trees/88_otus_img_gg_18may2012.tre
	filter_tree.py -i img_gg_otus_18may2012/trees/88_otus_img_gg_18may2012.tre -f img_gg_otus_18may2012/rep_set/85_otus_img_gg_18may2012.fasta -o img_gg_otus_18may2012/trees/85_otus_img_gg_18may2012.tre
	filter_tree.py -i img_gg_otus_18may2012/trees/85_otus_img_gg_18may2012.tre -f img_gg_otus_18may2012/rep_set/82_otus_img_gg_18may2012.fasta -o img_gg_otus_18may2012/trees/82_otus_img_gg_18may2012.tre
	filter_tree.py -i img_gg_otus_18may2012/trees/82_otus_img_gg_18may2012.tre -f img_gg_otus_18may2012/rep_set/79_otus_img_gg_18may2012.fasta -o img_gg_otus_18may2012/trees/79_otus_img_gg_18may2012.tre
	filter_tree.py -i img_gg_otus_18may2012/trees/79_otus_img_gg_18may2012.tre -f img_gg_otus_18may2012/rep_set/76_otus_img_gg_18may2012.fasta -o img_gg_otus_18may2012/trees/76_otus_img_gg_18may2012.tre
	filter_tree.py -i img_gg_otus_18may2012/trees/76_otus_img_gg_18may2012.tre -f img_gg_otus_18may2012/rep_set/73_otus_img_gg_18may2012.fasta -o img_gg_otus_18may2012/trees/73_otus_img_gg_18may2012.tre
	filter_tree.py -i img_gg_otus_18may2012/trees/73_otus_img_gg_18may2012.tre -f img_gg_otus_18may2012/rep_set/70_otus_img_gg_18may2012.fasta -o img_gg_otus_18may2012/trees/70_otus_img_gg_18may2012.tre


Acknowledgements
----------------

We wish to thank Amazon Web Services for the AWS in Education researcher's grant to the QIIME development group which supported to compute time used to develop this workflow and develop these reference collections.
