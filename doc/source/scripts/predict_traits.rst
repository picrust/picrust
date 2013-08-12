.. _predict_traits:

.. index:: predict_traits.py

*predict_traits.py* -- Given a tree and a set of known character states (observed traits and reconstructions), output predictions for unobserved character states
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


This script produces predictions of unobserved traits given a phylogenetic tree and a table that summarizes which traits are present in which ancestral organisms.
In the most common usage, this script is used to predict which gene families are present in each OTU (Operational Taxonomic Unit; roughly equivalent to a bacterial 'species'), given a tree and a set of ancestral state reconstructions.

The output of the script is a trait prediction file, which summarizes the predicted traits of each organism of interest (by default, this is all of the organisms that are tips in the phylogenetic tree).

The prediction method works as follows:

    1.  For each terminal (tip) node where a prediction is to be performed, the algorithm through the reconstructed ancestral states, and finds the last node in the ancestry of our organism of interest for which a prediction is available

    2.  The trait for the organism is then predicted based on a branch-length weighted average of the ancestral node and it's close relatives. (This is necessary because technical limitations involving the handling of ambiguous characters in many Ancestral State Reconstruction programs prevent the parent node of the organism from being directly reconstructed in most cases.)

    The exact weight function to use can be specified from the commandline (see options below).

    In general, this approach causes the prediction to be a weighted average of the closest reconstructed ancestor, and the either reconstructed or directly observed trait value of the organism of interest's sibling node(s).   



**Usage:** :file:`predict_traits.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-observed_trait_table
		The input trait table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format
	-t, `-`-tree
		The full reference tree, in Newick format
	
	**[OPTIONAL]**
		
	-o, `-`-output_trait_table
		The output filepath for trait predictions [default: predicted_states.tsv]
	-a, `-`-calculate_accuracy_metrics
		If specified, calculate accuracy metrics (i.e. how accurate does PICRUSt expect its predictions to be?) and add to output file [default: False]
	`-`-output_accuracy_metrics_only
		If specified, calculate accuracy metrics (e.g. NSTI), output them to this filepath, and do not do anything else. [default: None]
	-m, `-`-prediction_method
		Specify prediction method to use.  The recommended prediction method is set as default, so other options are primarily useful for control experiments and methods validation, not typical use.  Valid choices are:asr_and_weighting,nearest_neighbor,asr_only,weighting_only,random_neighbor.  "asr_and_weighting"(recommended): use ancestral state reconstructions plus local weighting with known tip nodes.  "nearest_neighbor": predict the closest tip on the tree with trait information.  "random_annotated_neighbor": predict a random tip on the tree with trait information. "asr_only": predict the traits of the last reconstructed ancestor, without weighting. "weighting_only": weight all genomes by distance, to the organism of interest using the specified weighting function and predict the weighted average.   [default: asr_and_weighting]
	-w, `-`-weighting_method
		Specify prediction the weighting function to use.  This only applies to prediction methods that incorporate local weighting ("asr_and_weighting" or "weighting_only")  The recommended weighting  method is set as default, so other options are primarily useful for control experiments and methods validation, not typical use.  Valid choices are:exponential,linear,equal.  "exponential"(recommended): weight genomes as a negative exponent of distance.  That is 2^-d, where d is the tip-to-tip distance from the genome to the tip.  "linear": weight tips as a linear function of weight, normalized to the maximum possible distance (max_d -d)/d. "equal_weights": set all weights to a constant (ignoring branch length).   [default: exponential]
	-l, `-`-limit_predictions_by_otu_table
		Specify a valid path to a legacy QIIME OTU table to perform predictions only for tips that are listed in the OTU table (regardless of abundance)
	-g, `-`-limit_predictions_to_organisms
		Limit predictions to specific, comma-separated organims ids. (Generally only useful for lists of < 10 organism ids, for example when performing leave-one-out cross-validation).
	-r, `-`-reconstructed_trait_table
		The input trait table describing reconstructed traits (from `ancestral_state_reconstruction.py <./ancestral_state_reconstruction.html>`_) in tab-delimited format [default: None]
	`-`-confidence_format
		The format for the confidence intervals from ancestral state reconstruction. Only needed if passing a reconstruction confidence file with -c or --reconstruction_confidence.  These are typically sigma values for maximum likelihood ASR  methods, but 95% confidence intervals for phylogenetic independent contrasts (e.g. from the ape R packages ace function with pic as the reconstruction method).  Valid choices are:sigma,confidence_interval. [default: sigma]
	-c, `-`-reconstruction_confidence
		The input trait table describing confidence intervals for reconstructed traits (from `ancestral_state_reconstruction.py <./ancestral_state_reconstruction.html>`_) in tab-delimited format [default: None]
	`-`-output_precalc_file_in_biom
		Instead of outputting the precalculated file in tab-delimited format (with otu ids as row ids and traits as columns) output the data in biom format (with otu as SampleIds and traits as ObservationIds) [default: False]


**Output:**

Output is a table (tab-delimited or .biom) of predicted character states


Required options with NSTI:

::

	predict_traits.py -a -i trait_table.tab -t reference_tree.newick -r asr_counts.tab -o predict_traits.tab

Limit predictions to particular tips in OTU table:

::

	predict_traits.py -a -i trait_table.tab -t reference_tree.newick -r asr_counts.tab -o predict_traits_limited.tab -l otu_table.tab

Reconstruct confidence

::

	predict_traits.py -a -i trait_table.tab -t reference_tree.newick -r asr_counts.tab -c asr_ci.tab -o predict_traits.tab


