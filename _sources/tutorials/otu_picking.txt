.. _otu_picking_tutorial:

Picking OTUs for use in PICRUSt
===============================

Introduction
------------
This document covers how to pick OTUs from marker gene data to use with PICRUSt. To do this, you'll use a 'closed-reference' OTU picking protocol where you search sequences against the GG reference OTUs at a specified percent identity, and discard any reads that don't hit that reference collection. The newest available reference collection can be found here:

 * gg_13_5_otus.tar.gz (`download <http://greengenes.secondgenome.com/downloads/database/13_5>`_ )

This tutorial assumes that you have QIIME installed. See the `QIIME website <http://www.qiime.org>`_ for details on how to use QIIME. The quickest way to get started with QIIME is working on the Amazon Web Services cloud, and you can find instructions for `using QIIME on the cloud here <http://qiime.org/tutorials/working_with_aws.html>`_.

Picking closed reference OTUs with QIIME
----------------------------------------

To pick 'closed reference' OTUs with QIIME for use in PICRUST, you should begin with a demuliplexed fasta file in QIIME format, and the IMG/GG reference collection.

The demuliplexed fasta file in QIIME format is a standard multi-line fasta file, where sequence identifiers are of the form ``sampleID_seqID``. In these sequence identifiers, ``sampleID`` indicates the name of the sample, and ``seqID`` is a unique (with respect to the current file) sequence identifier. The ``seqID`` values are often just assigned as ascending integers to ensure their uniqueness. These files can be generated with the ``split_libraries`` * scripts in QIIME, or using other software. This file might look like::

	>s1_1
	ACCTTAGGATGGATTAGACCCAGA
	>s1_2
	ACCAGGATGGACCCCTTAGACCCAGA
	>s2_3
	AGGTCCTGGATGGACCCCTTAGACCCAGA
	>s1_4
	ACCAGTTGATGGACCCCTTAGACCCAGA

In this example we have two samples (``s1`` and ``s2``) and four sequences, three of which are in ``s1`` and one of which is in ``s2``. For additional details on this file, see `here <http://qiime.org/documentation/file_formats.html#demultiplexed-sequences>`_.

If your demultiplexed fna file is named ``seqs.fna`` and is in your current working directory, and you have unzipped the GG OTUs in the current working directory, you would pick OTUs with the following commands::

	echo "pick_otus:enable_rev_strand_match True"  >> $PWD/otu_picking_params_97.txt
	echo "pick_otus:similarity 0.97" >> $PWD/otu_picking_params_97.txt
	pick_closed_reference_otus.py -i $PWD/seqs.fna -o $PWD/ucrC97/ -p $PWD/otu_picking_params_97.txt -r $PWD/gg_13_5_otus/rep_set/97_otus.fasta

This command picks OTUs and builds a `biom-formatted OTU table <http://www.biom-format.org>`_, with OTUs assigned at 97% identity. The primary file of interest will be ``ucrC97/uclust_ref_picked_otus/otu_table.biom``, which will be the OTU table that you pass to PICRUSt. 

Alternatively, if you'd like to pick OTUs at a 90% percent identity, you could run the following command::

	echo "pick_otus:enable_rev_strand_match True"  >> otu_picking_params_90.txt
	echo "pick_otus:similarity 0.90" >> otu_picking_params_90.txt
	pick_closed_reference_otus.py -i $PWD/seqs.fna -o $PWD/ucrC90/ -p $PWD/otu_picking_params_90.txt -r $PWD/gg_13_5_otus/rep_set/97_otus.fasta

If you have parallel QIIME running in your environment, you can modify the commands to run in parallel by append the ``-a`` flag and the ``-O`` parameter. This might looking like the following to use 8 processors to pick OTUs at 97% identity::

	echo "pick_otus:enable_rev_strand_match True"  >> $PWD/otu_picking_params_97.txt
	echo "pick_otus:similarity 0.97" >> $PWD/otu_picking_params_97.txt
	pick_closed_reference_otus.py -i $PWD/seqs.fna -o $PWD/ucrC97/ -p $PWD/otu_picking_params_97.txt -r $PWD/gg_13_5_otus/rep_set/97_otus.fasta -a -O 8

Note that you should not specify a value for ``-O`` that is greater than the number of processor/cores that you have available. This will result in longer run times as multiple jobs will compete for the same processor.

