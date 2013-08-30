.. _qiime_tutorial:

Analyzing metagenomes with QIIME
================================

.. include:: ../global.rst

Because the metagenomes are provided in BIOM format by default, these can be plugged into many of the downstream analysis tools available in `QIIME <www.qiime.org>`_. QIIME's `Shotgun Metagenome Analysis tutorial <http://qiime.org/tutorials/shotgun_analysis.html>`_ illustrates a couple of the steps that can be applied. The steps that will primarily be of interest in that tutorial are the ones that begin with a ``.biom`` file. For example, `computing beta diversity and PCoA plots <http://qiime.org/tutorials/shotgun_analysis.html#computing-beta-diversity-and-generating-pcoa-plots>`_ and `generating summaries of samples by KO categories <http://qiime.org/tutorials/shotgun_analysis.html#generating-summaries-of-samples-by-ko-category>`_. 

Many of `QIIME's tutorials that describe diversity analyses <http://qiime.org/tutorials/index.html>`_ are applicable to PICRUSt-predicted metagenome tables. Specific analysis tools that may be useful include:

	* `alpha_diversity.py <http://qiime.org/scripts/alpha_diversity.html>`_
	* `beta_diversity.py <http://qiime.org/scripts/beta_diversity.html>`_
	* `compute_core_microbiome.py <http://qiime.org/scripts/compute_core_microbiome.html>`_
	* `jackknifed_beta_diversity.py <http://qiime.org/scripts/jackknifed_beta_diversity.html>`_
	* `make_distance_boxplots.py <http://qiime.org/scripts/make_distance_boxplots.html>`_
	* `alpha_rarefaction.py <http://qiime.org/scripts/alpha_rarefaction.html>`_
	* `beta_diversity_through_plots.py <http://qiime.org/scripts/beta_diversity_through_plots.html>`_
	* `otu_category_significance.py <http://qiime.org/scripts/otu_category_significance.html>`_
	* `shared_phylotypes.py <http://qiime.org/scripts/shared_phylotypes.html>`_


Plots of functional categories at various levels can be created using `summarize_taxa_through_plots.py <http://qiime.org/scripts/summarize_taxa_through_plots.html>`_
	
	* Since KEGG Orthologs belong to several pathways you should collapse your PICRUSt predictions to the desired hierarchy level using :ref:`categorize_by_function.py <categorize_by_function>` ::

	        categorize_by_function.py -i metagenome_predictions.biom -c "KEGG_Pathways" -l 2 -o metagenome_at_level2.biom	
	
	* Then add the following lines to a `qiime parameter file <http://qiime.org/documentation/qiime_parameters_files.html>`_ (e.g. qiime_params.txt) ensuring that the level you collapsed at is the same in your config file ::
	       
	        summarize_taxa:md_identifier    "KEGG_Pathways"
		summarize_taxa:absolute_abundance   True 
		summarize_taxa:level    2

	* Lastly, run `summarize_taxa_through_plots.py <http://qiime.org/scripts/summarize_taxa_through_plots.html>`_ ::

	        summarize_taxa_through_plots.py -i metagenome_at_level2.biom -p qiime_params.txt -o plots_at_level2


There are also a number of scripts in QIIME that may be useful for more general processing of your BIOM table. These include the following:

	* `single_rarefaction.py <http://qiime.org/scripts/single_rarefaction.html>`_
	* `filter_otus_from_otu_table.py <http://qiime.org/scripts/filter_otus_from_otu_table.html>`_
	* `filter_samples_from_otu_table.py <http://qiime.org/scripts/filter_samples_from_otu_table.html>`_
	* `per_library_stats.py <http://qiime.org/scripts/per_library_stats.html>`_
	* `filter_taxa_from_otu_table.py <http://qiime.org/scripts/filter_taxa_from_otu_table.html>`_
	* `merge_otu_tables.py <http://qiime.org/scripts/merge_otu_tables.html>`_
	* `sort_otu_table.py <http://qiime.org/scripts/sort_otu_table.html>`_
	* `split_otu_table.py <http://qiime.org/scripts/split_otu_table.html>`_
	* `split_otu_table_by_taxonomy.py <http://qiime.org/scripts/split_otu_table_by_taxonomy.html>`_

Note that while many of these refer to OTU table, it's just a nomenclature issue. These are generally applicable to ``.biom`` tables.

Finally, if you're interested in comparing real to predicted metagenomes, or predicted metagenomes to 16S, you'll be interested in the `Procrustes Analysis tutorial <http://qiime.org/tutorials/procrustes_analysis.html>`_ and the `Comparing Distance Matrices tutorial <http://qiime.org/tutorials/distance_matrix_comparison.html>`_.
