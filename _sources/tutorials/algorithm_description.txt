.. _algorithm_description:

How PICRUSt Works
=================

.. include:: ../global.rst

Overview
---------
PICRUSt is designed to estimate the gene families contributed to a metagenome by bacteria or archaea identified using 16S rRNA sequencing.  Intermediate steps in this pipeline may also be of independent interest, as they allow for phylogenetic prediction of organismal traits using reference examples (here applied to the problem of gene content prediction), and correction for variable marker gene copy number.  

For an alternative presentation of the same algorithm, along with controls please see the `PICRUSt paper <http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.2676.html>`_


Gene Content Prediction
-----------------------

A common and heavily studied problem in phylogenetics involves estimating the properties of ancestral organisms from living relatives (ancestral state reconstruction or ASR).  PICRUSt extends these methods to predict traits of living organisms that are unknown not because they belong to ancestral organisms, but because they have simply not yet been studied.   The main application of this Hidden State Prediction (HSP) in the core PICRUSt workflow is to estimate the gene content of microorganisms for which no genome sequence is available, by using their sequenced relatives as a reference. 

The central idea is that despite various important forms of microbial genome plasticity (gene loss, duplication, or gene transfer), the genes present in microbial genomes are much more similar amongst related bacteria or archaea than distant relatives.  Therefore, when sufficient genome sequences are available, it is possible to predict which gene families are present  in a given microbial OTU from phylogeny alone.  The important caveat is the 'when sufficient genome sequences are available'.   Because of the capability of microbial genomes to change rapidly over evolutionary time, there will always be uncertainty in these predictions.   OTUs lacking closely related genomes will be less easily predicted, as will gene families that undergo gene loss or duplication rapidly over evolutionary time (for example,  because they are prone to HGT).  Therefore, each prediction is optionally accompanied by a 95% confidence interval that reflects these sources of uncertainty.  

*These gene content predictions are precalculated for protein-coding genes present in KEGG or COG gene families and 16S rRNA gene copy number.  Therefore users do not typically need to recalculate the gene content predictions to use PICRUSt*


Gene content prediction takes two input data files:  a table of traits to be predicted, and a phylogenetic tree (see the :ref:`Gene Content Prediction tutorial <genome_prediction_tutorial>` for details on the file format).  Typically, the table of traits is a table of gene copy numbers for each gene family in each sequenced bacterial and archaeal genome in the IMG database, and the phylogenetic tree is the Greengenes phylogeny.   However, alternative trees and trait tables could be employed to e.g. attempt to predict the traits of an animal based on those of phylogenetic relatives. 

1. Tree Pruning and Formatting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    **Rationale**: Ancestral state reconstruction (ASR) software generally requires that ids match exactly between the matrix of character states (i.e. organismal traits) and the phylogeny.   Therefore PICRUSt must generate a `pruned` version of the full tree before the input data can be provided to ASR methods. 
    
    **Details**: PICRUSt's format_tree_and_trait_table.py script generates a pruned tree and reduced trait table that contain only organims shared between the tree and the trait table.   This step also formats trees for input into ASR programs.   This includes:  forcing trees to be fully bifurcating (polytomies are resolved arbitrarily with epsilon length branches), substituting any zero-length branches with a small epsilon branch length, and modifying some ids to avoid characters that would break ASR software.  See the format_tree_and_trait_table.py documentation for all available options.The phylogeny provided is currently assumed to be correct.  Therefore, additional investigations of the effects of error during phylogenetic reconstruction, if desired, must be conducted by using alternative inputs to PICRUST (e.g. of trees derived from bootstrapped versions of the phylogeny).      

2. Ancestral State Reconstruction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    **Rationale**: Because much gene content is conserved between closely related microorganisms, an accurate picture of the genes present in the common ancestor of an organism and a sequenced genome can provide a baseline prediction for the gene content in a modern organism that we have not yet studied.  Therefore, PICRUSt uses one of several ancestral state reconstruction algorithms to estimate ancestral gene content.
    
    **Details**: PICRUSt's ancestral_state_reconstruction.py script applies one of several external ancestral state reconstruction programs to generate ancestral states.   When maximum likelihood methods are employed, the count of each gene family is treated as a continuous evolutionary character evolving under a Brownian Motion model.  By default, a fast approximate ASR method (ace_pic from the ape R package) but a full maximum likelihood inference is available as a commandline option.   Treating gene copy number as a continuous character bears some explaination.  An alternative approaches might treat gene counts as discrete states.  This perhaps seems more natural, and has the potential benefit of allowing asymmetrical rate matrices.  However, gene copy numbers for some gene families can reach into the hundreds in certain genomes.  These high copy numbers produce extremely large rate matrices that may both render ancestral state reconstruction computationally intractable, and may also not be sufficiently contrained by available data.  Therefore the simpler continuous model of evolution is currently employed.

    Ancestral state reconstruction under a Brownian motion, continuous-trait model also produces an estimate of the uncertainty in each ancestral state. Finally, one rate (Brownian motion parameter, Sigma) is produced for each character (e.g. gene family).  This rate reflects how rapidly that gene family changes in copy number over evolutionary time. This rate therefore implicitly incorporates uncertainty due to a variety of mechanisms of microbial genome plasticity, including horizontal gene transfer.  These uncertainty estimates and rates of gene evolution are optionally used during trait prediction to propagate uncertainty in the ancestral state reconstruction forward to the gene prediction. 


3. Predict Traits
^^^^^^^^^^^^^^^^^  
   **Rationale**:  Following ancestral state reconstruction, PICRUSt shoud have available three key pieces of information: the traits of studied organisms (e.g. gene counts in each sequenced genome), the traits of ancestral organisms, and  a reference tree containing both the reference organisms and the inferred ancestral organisms.  The estimate for the gene content of an organism lacking a genome sequence  will depend on its position in the phylogenetic tree relative to organims with sequenced genomes and/or inferred ancestral genomes.   PICRUST estimates the  value of these traits by weighting sequenced or ancestral genomes based on their distance.  This weighting can be linear, exponential, or based on variance.  Alternatively, trait prediction can optionally be accomplished without reference to ancestral states by simply predicting the traits of the most closely related genome on the supplied phylogeny.

   **Details**: Optionally, trait prediction can output a confidence interval for each predicted trait. This confidence interval is based on the variance of the prediction.  This variance is based on two factors:  the uncertainty in **(1)** each ancestral state reconstruction itself, and **(2)** the additional uncertainty due to evolutionary distance separating the organism for which traits are being predicted from the reference organisms and/or reconstructed ancestors that drive the prediction.  
   
   The first type of variance is accounted for during trait prediction by standard formulas for propagating variance across a weighted average (see the `predict_traits.py <predict_traits>` library itself for full details). Variance of the second type, owing to recent evolutionary changes, is accounted for by using the rate parameter inferred during ancestral state reconstruction. This has the effect that long branches (relative to the LCA with a sequenced genome or ancestror for which a reconstructed state is available) or highly variable genes can both in theory result in high variances.  In practice, however, it appears that effects due to the availability or unavailability of reference genomes have a much greater effect than fast- vs. slow- evolving gene families. These results and some benchmarks for the 95% confidence intervals are discussed in the Supplemental material to the `PICRUSt paper <http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.2676.html>`_


16S rRNA Copy Number Prediction
-------------------------------
    
    **Rationale**:The 16S rRNA gene copy number is predicted in the same fashion as other genes.  However, because counts of 16S rRNA gene operons per genome are derived from different source data in IMG, prediction for 16S rRNA copy number is run separately.
    
    **Details**: In order to account for variable 16S rRNA operon copy numbers between bacterial taxa, 16S rRNA copy numbers are predicted by predict_traits.py in exactly the same manner as other genes, but using the IMG counts of 16S rRNA copy number per genome as the input.  These predicted copy numbers are available as a precalculated file :ref:`picrust_precalculated_files`.  Benchmarks for 16S operon copy number prediction accuracy can be found in the Supplemental materials to the `PICRUSt paper <http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.2676.html>`_.  

Copy Number Normalization
-------------------------
    
    **Rationale**: Because 16S rRNA gene operons can vary from 1 to 15 copies in bacteria, observed relative abundances in 16S rRNA sequencing studies may vary from true organismal abundances.  For example, the metagenomic contribution of a bacterium with 4 rRNA operon copies might be inflated by 4 fold if no copy number normalization was applied.  For further discussion of 16S copy number normalization, see `Kembel et al 2012 <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002743>`_

    **Details**: :ref:`normalize_by_copy_number.py <normalize_by_copy_number>` simply divides each OTU value in each sample by its predicted 16S rRNA copy number in the :ref:`picrust_precalculated_files`.   Currently, variance in the metagenomic confidence interval is not adjusted due to uncertainty in 16S rRNA copy number estimation, but this is a straightforward addition for future development.  Benchmarks for 16S operon copy number prediction accuracy can be found in the Supplemental materials to the `PICRUSt paper <http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.2676.html>`_.  Although the accuracy of 16S rRNA predictions has has been tested, and there are clear reasons to think this might affect overall accuracy, the overall effect of copy number inference on beta-diversity clustering or PICRUSt metagenome prediction is an area that would benefit from further characterization.



Metagenome Prediction
---------------------
    
    **Rationale**:   Following copy number normalization, the OTU table should now correspond to the relative abundances of organisms rather than the relative abundance of the 16S rRNA itself.  Because the number of gene copies for each gene family per organism has already been estimated in the precalculated files produced by :ref:`predict_traits.py <predict_traits>`, producing a metagenome prediction is handled by simply multiplying the vector of gene counts for each OTU by the abundance of that OTU in each each sample, and summing across all OTUs.  
    
    **Details**:  For computational efficiency, these sums are handled using numpy's dot product. Optionally, uncertainty in trait prediction  (see above) can be propagated through trait prediction to produce a 95% confidence interval for the prediction. This is done by breaking down the individual steps involved in calculating the metagenome prediction, and propagating variance across each. Unlike the confidence intervals for specific per-genome predictions, PICRUSt's capability for 95% confidence intervals on metagenome predictions, although unit-tested, has not yet been validated on biological datasets.  Therefore metagenomic confidence intervals should currently be considered experimental.   See :ref:`predict_metagenomes.py <predict_metagenomes>` and associated library code for full details.    
