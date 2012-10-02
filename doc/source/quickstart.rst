.. _quickstart:

PICRUST Quickstart Guide
========================

This gives the barebone commands needed to run PICRUST from start to finish.

.. include:: global.rst

1. Download PICRUST Software
----------------------------
* Download the `PICRUST development software`_ using SVN:: 
	
	svn checkout svn://svn.code.sf.net/p/picrust/code/trunk ~/picrust-dev

* Install PICRUST::

	#For Bash shell 
	echo 'export PYTHONPATH=~/picrust-dev/picrust/:$PYTHONPATH' | cat >> ~/.bashrc	
	echo 'export PATH=~/picrust-dev/scripts/:$PATH' | cat >> ~/.bashrc	

	OR

	#For tsch shell
	echo 'setenv PYTHONPATH ~/picrust-dev/picrust/:$PTYHONPATH' | cat >> ~/.tschrc
	echo 'setenv PATH ~/picrust-dev/scripts/:$PATH' | cat >> ~/.tschrc

2. Download PICRUST Data Files
------------------------------

* Download the `PICRUST precalculated files`_::

	cd ~/picrust-dev && wget -O - "http://dl.dropbox.com/s/f8sy84kxi2havww/picrust_precalculated_files.tgz" | tar -xvf -
	

3. Prepare your OTU table
-------------------------

* Download the `PICRUST GG reference data`_::
	
	cd ~/picrust-dev && wget -O - "http://s3.amazonaws.com/picrust-public-data/img_gg_otus_18may2012.tgz" | tar -xzf -

* Pick reference based OTUs (assuming QIIME already installed and demultiplexed fasta file (seqs.fna))::

	cat "pick_otus:enable_rev_strand_match True"  >> $PWD/otu_picking_params_97.txt
	cat "pick_otus:similarity 0.97" >> $PWD/otu_picking_params_97.txt
	pick_reference_otus_through_otu_table.py -i $PWD/seqs.fna -o $PWD/ucrC97/ -p $PWD/otu_picking_params_97.txt -r $PWD/img_gg_otus_18may2012/rep_set/97_otus_img_gg_18may2012.fasta

* For more information see :ref:`otu_picking_tutorial`.


4. Run PICRUST
--------------

* Normalize your OTU table::

	normalize_by_copy_number.py -i $PWD/ucrC97/uclust_ref_picked_otus/otu_table.biom -c ~/picrust-dev/picrust_precalculated_files/16S_acepic_predict_traits_97.biom.gz -o your_normalized_otu_table.biom

* Get KEGG predictions (NOTE: This step currently requires approx 5GB RAM)::

	predict_metagenomes.py -i your_normalized_otu_table.biom -c ~/picrust-dev/picrust_precalculated_files/KEGG_acepic_predict_traits_97.biom.gz -o your_KEGG_predictions.biom

* For more information see :ref:`metagenome_prediction_tutorial`
