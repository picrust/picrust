.. include:: ../global.rst
.. _quickstart_guide:

Quickstart Guide
================

This gives a basic overview of how to use PICRUSt. Before starting, you should have installed PICRUSt as described in :ref:`install`.

Prepare your OTU table
-------------------------

* Download the `PICRUSt GG reference data`_::
	
	cd ~/picrust-dev && wget -O - "http://s3.amazonaws.com/picrust-public-data/img_gg_otus_18may2012.tgz" | tar -xzf -

* Pick reference based OTUs (assuming QIIME already installed and demultiplexed fasta file (seqs.fna))::

	echo "pick_otus:enable_rev_strand_match True"  >> $PWD/otu_picking_params_97.txt
	echo "pick_otus:similarity 0.97" >> $PWD/otu_picking_params_97.txt
	pick_reference_otus_through_otu_table.py -i $PWD/seqs.fna -o $PWD/ucrC97/ -p $PWD/otu_picking_params_97.txt -r $PWD/img_gg_otus_18may2012/rep_set/97_otus_img_gg_18may2012.fasta

* For more information see :ref:`otu_picking_tutorial`.


Run PICRUSt
--------------

* Normalize your OTU table (:ref:`normalize_by_copy_number`)::

	normalize_by_copy_number.py -i $PWD/ucrC97/uclust_ref_picked_otus/otu_table.biom -o normalized_otus.biom

* Get KEGG predictions (:ref:`predict_metagenomes`) (**NOTE: This step currently requires approx 5GB RAM**)::

	predict_metagenomes.py -i normalized_otus.biom -o metagenome_predictions.biom

* For more information see :ref:`metagenome_prediction_tutorial`
