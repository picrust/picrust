.. _genome_prediction_tutorial:

Inferring Individual Genomes
============================

Introduction
------------

This tutorial explains the initial steps in the **PICRUSt** pipeline that looks after formatting input tree and trait tables, running ancestral state reconstruction, and extending predictions to all tips in the reference tree.  

Essential Files
---------------
All files required run this tutorial are here (https://www.dropbox.com/s/l9djlyssn3cspz8/genome_prediction_tutorial_files.zip). Descriptions of these files are below. 

* Reference Tree (.nwk)
    * This is a phylogenetic tree constructed from a marker gene (typically 16S), that has tips representing both sequenced genomes and non-sequenced genomes. 
    * We will use the greengenes tree **gg_tree.nwk**.

* Marker gene copy number (tab-delimited .tab)
    * This is a simple tab-delimited table with two columns. The first column contains the genome identfiers column ids (typically IMG), and the second column contains the number of copies of the marker gene (16S). The first row must contain header ids for the columns.
    * We will use the file **IMG_16S_counts.tab**.

* Functional trait copy number (tab-delimited .tab)
    * Same format as the marker gene copy number table above with genome identifiers being the first column, but contains an additional column for every functional trait with the column id representing that trait (e.g. K00001).
    * We will use KEGG KO functions **IMG_ko_counts.tab**.

* Mapping file for tree tip ids to genome ids (tab-delimited .tab)
    * Simple two column tab delimited file with first column containing ids of tips in the tree that have genome sequences. Second column contains genome identifiers from marker gene and functional trait copy number files. 
    * We will use the file **GG_to_IMGv350.txt**.
    * *(Note: This file is optional if the tips representing genomes in the reference tree exactly match genome identifiers in the trait tables)*. 


Formatting tree and trait tables
--------------------------------
`format_tree_and_trait_table.py <../scripts/format_tree_and_trait_table.html>`_ does numerous formatting and checks to the reference tree and the trait tables. 
The following steps are done for both the marker gene copy number table and functional trait copy number table. 

1. All internal nodes in the reference tree are checked for problematic characters and unlabelled internal nodes are given labels. 
2. A pruned tree is created that contains only tips that have copy number predictions from sequenced genomes.
3. Any traits in the trait table that are not in the reference tree are removed. 

Format the 16S marker gene copy number table: ::

	format_tree_and_trait_table.py -t GG_tree.nwk -i IMG_16S_counts.tab -m GG_to_IMGv350.txt -o format/16S/

Format the KO functional trait copy number table: ::

	format_tree_and_trait_table.py -t GG_tree.nwk -i IMG_ko_counts.tab -m GG_to_IMGv350.txt -o format/KEGG/

Each of the above two commands creates 3 output files in the directory specified by the '-o' option. 

1. **reference_tree.newick**
2. **pruned_tree.newick**
3. **trait_table.tab**

Note: that the reference_tree.newick files from the two commands will be identical. 

Ancestral State Reconstruction
------------------------------
`ancestral_state_reconstruction.py <../scripts/ancestral_state_reconstruction.html>`_ runs ancestral state reconstruction (ASR) to make predictions for each trait for every internal node in the pruned tree. 

Various methods of ASR can be chosen including 'wagner', 'ace_pic', 'ace_ml', and 'ace_reml'. The first two methods are very fast while the latter two ML methods require parallelization to be finished in a suitable timeframe. 

Input is a formatted trait table and a pruned tree.

First, we will use the default ASR method 'wagner' on the 16S trait table and the pruned tree: ::

	ancestral_state_reconstruction.py -i format/16S/trait_table.tab -t format/16S/pruned_tree.newick -o asr/16S_asr_wagner_counts.tab 

Second, run the same method on the formatted KO trait table (with the pruned tree): ::

	ancestral_state_reconstruction.py -i format/KEGG/trait_table.tab -t format/KEGG/pruned_tree.newick -o asr/KEGG_asr_wagner_counts.tab

Output is a tab delimited file with first column containing internal node labels from the tree and other columns containing predicted ASR counts for each trait.

(Optional) Running other ASR methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Other ASR methods can be run using the '-m' option: ::

	ancestral_state_reconstruction.py -i format/16S/trait_table.tab -t format/16S/pruned_tree.newick -o asr/16S_asr_acepic_counts.tab -m ace_pic

ML based ASR methods (`ace_reml` and `ace_ml`) require much longer computational times. Parallelization is supported for multiple platforms (e.g. SGE, Torque, etc.). 

For example use '-p' option to turn on parallelization, '-j sge' to specify a SGE based cluster, and '-n 200' to group commands into 200 jobs: ::

	ancestral_state_reconstruction.py -i format/KEGG/trait_table.tab -t format/KEGG/pruned_tree.newick -o asr/KEGG_asr_aceml_counts.tab -m ace_ml -p -j sge -n 200

Genome (trait) Prediction
-------------------------
`predict_traits.py <../scripts/predict_traits.html>`_ extends ASR predictions from internal nodes and real values from known genomes to tips in the genome tree with no genome information.

Input is the formatted trait table, the formatted reference tree, and the ASR output file.

Make predictions for 16S: ::

	predict_traits.py -i format/16S/trait_table.tab -t format/16S/reference_tree.newick -r asr/16S_asr_wagner_counts.tab -o predict_traits/trait_predictions_16S_wagner.biom 

Make predictions for KOs: ::
	
	predict_traits.py -i format/KEGG/trait_table.tab -t format/KEGG/reference_tree.newick -r asr/KEGG_asr_wagner_counts.tab -o predict_traits/trait_predictions_KEGG_wagner.biom


Output is a biom formatted file with 'Observations' (like rows) as tree tip ids (e.g. genomes/OTUs) and 'Samples' (like columns) as functional traits. 

(Optional) Limiting predictions to those in OTU table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`predict_traits.py` can take a long time to run if making predictions for all tips in the green genes reference tree (400k tips). Therefore, you can limit the number of predictions to only those in your metagenome OTU table (the ones you care about) using the '-l' option. 

Make predictions for KOs for a given OTU table using '-l' option: ::
	
	predict_traits.py -i format/KEGG/trait_table.tab -t format/KEGG/reference_tree.newick -r asr/KEGG_asr_wagner_counts.tab -l your_otu_table.tsv -o predict_traits/your_otu_trait_predictions_KEGG_wagner.biom 


