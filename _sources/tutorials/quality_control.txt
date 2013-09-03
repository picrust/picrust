.. _quality_control:

Quality Control of PICRUST Predictions
======================================

.. include:: ../global.rst

This section covers steps that can be taken to ensure PICRUSt predictions will be as accurate as possible, and to characterize how well or poorly the gene content of a given set of 16S rRNA samples can be predicted.

Overview
--------

PICRUSt’s main use is in estimating the bacterial and archaeal genes present in a microbial community metagenome, using 16S rRNA data.   The output is a table of gene family abundances for each sample, with the exact scheme used to delineate gene families chosen by you (options include KEGG, COG, PFAM and RFAM).  The accuracy of PICRUSt’s predicted metagenomes can vary quite a bit depending on several factors including, most importantly, the extent to which sequenced genomes are available for the most abundant species in your community of interest.   

More detail about the types of tests we have done on the algorithm, and the detailed results are available in the `PICRUSt manuscript in Nature Biotechnology <http://www.nature.com/nbt/journal/vaop/ncurrent/abs/nbt.2676.html>`_

Here, we will cover the most important conclusions from that research, and delve more into some extra steps you can take in order to predict how accurate PICRUSt will be for your samples of interest.  First we will discuss some :ref:`general limitations <general_limitations>` of PICRUSt that may be worth considering when interpreting results.  Second, we present quality control steps that can be taken to quantify or predict PICRUSt's accuracy and determine if it will be useful for your dataset:

The :ref:`basic quality control checklist <basic_qc_checklist>` include tests that you can run along with your PICRUSt analysis without doing any additional sequencing.  

The :ref:`empirical quality control checklist <empirical_qc_checklist>` is intended for larger projects (with for example hundreds or thousands of 16S rRNA libraries and perhaps a dozen or more shotgun metagenomes available), and includes validations similar to those used in the PICRUSt manuscript.   At the bottom of the page, we include a checklist of items you may wish to include in your methods section when describing PICRUSt so that others can duplicate your findings. 

Both the ‘basic’ and ‘advanced’ checklists assume you are doing the main PICRUSt workflow- predicting a ‘virtual metagenome’ from 16S rRNA data.  However, PICRUSt can in theory be used to predict any trait (continuous evolutionary character) across organisms.   Additional discussion of quality control steps you might take if you are using PICRUSt in a customized workflow will be discussed in a forthcoming ‘custom trait prediction with PICRUSt’ section.

Finally, we discuss some ideas for :ref:`reporting a PICRUSt analysis <reporting_picrust_results>` such that it can be reproduced by others.

We would like this page to be a resource for the community.  If you have additional suggestions for additional steps that will help the community with their PICRUSt analysis, we would love to hear them-  please e-mail picrust-users@googlegroups.com

.. _general_limitations:

General limitations of PICRUSt
------------------------------

* Since the input data for the standard PICRUSt workflow is 16S rRNA, any eukaryotic or viral contributions to the metagenome will not be predicted.  Therefore it is best to think of PICRUSt as predicting the *portion* of the full metagenome contributed by the organisms targeted by your primers.

* Biased primers may result in inaccurate predictions.  Only genes from organisms amplified by your primer will be included.  The main use for PICRUSt is in taking 16S rRNA data, and predicting a metagenome from that data using evolutionary modelling of how gene content has changed relative to sequenced genomes.   Therefore PICRUSt can only predict the portion of the metagenome that is contributed by the set of organisms picked up by your primers.   The PICRUSt validation datasets used universal 515f/806r V4 16S rRNA primers designed to minimize (though not eliminate) bias across bacterial/archaeal taxonomy [link Caporaso et al 2011, Walters et al 2011].   If the primers used don’t amplify an organism, then of course that organisms’ contribution to the metagenome is not predicted.   As an example, many popular 16S rRNA primers including 27F/338R did not work well for amplifying *Verrucomicrobia* (see `Bergmann et al 2011 <http://www.ncbi.nlm.nih.gov/pubmed/22267877>`_).   If your sample had a large proporition of *Verrucomicrobia* in it, and you used such a primer set, then of course the metagenome predicted by PICRUSt would also underestimate genes contributed by *Verrucomicrobia*.   
  
* PICRUSt can only predict gene families that are already known and included in the orthology reference used (KEGG KOs by default).  Therefore, genes that do not fall within the orthology scheme used will not be predicted, despite the possibility that such 'dark matter' could play an important role in the system at hand.    

* PICRUSt output that maps gene families to putative functions or pathways  is purely based on the particular input reference used. Therefore, any gaps or inaccuracies in pathway annotation or assignments of gene function will still be present.  As an example, many KEGG Orthology groups are listed as participating in pathways not found in bacteria or otherwise not reflective of true function.  In many cases this is simply due to bacteria containing  (distant) homologs of enzymes with important roles in, for example,  mammalian pathways.  Therefore, it is worth carefully checking KEGG pathway annotations to ensure that they are reasonable for your system. 

* PICRUSt is based on evolutionary modeling of the gene contents of known reference genomes.  Therefore, accuracy for any given sample type will depend heavily on the availability of appropriate references.  See the checklists below for some methods for quantifying whether this will be an issue for your samples, and for calculating 95% confidence intervals around gene content predictions.

.. _basic_qc_checklist:

Basic Quality Control steps for PICRUSt
---------------------------------------

These are some steps that can be taken to ensure that you get the most out of your PICRUSt analysis. 

1. Resources for understanding the PICRUSt algorithm and its performance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The PICRUSt algorithm and a number of accuracy tests are described in the `PICRUSt paper <http://www.nature.com/nbt/journal/vaop/ncurrent/abs/nbt.2676.html>`_. An additional description of the algorithm is available on the PICRUSt webpage here: :ref:`algorithm_description`.  Importantly, testing indicates that while PICRUSt can work quite well when adequete reference genomes are available, results may be unreliable for underexplored communities such as the Guerrero Negro microbial mats where adequete genome coverage is not currently available. 

2. Check input OTU tables
^^^^^^^^^^^^^^^^^^^^^^^^^

Appropriate quality control on the input 16S rRNA data is important for ensuring an accurate prediction when using PICRUSt.   For the standard metagenome prediction workflow, PICRUSt requires that input 16S rRNA data be in the form of a `BIOM format table <http://biom-format.org/index.html>`_  table that was picked by mapping your reads to references in the greengenes tree.   The most straightforward way to do this is using QIIME’s reference-based OTU picking pipeline as described in the :ref:`OTU Picking Tutorial <otu_picking_tutorial>`.   

 This has several implications:
        
* **OTU IDs must be based on greengenes**: If you input an OTU table generated some other way, but with some IDs that happen to match greengenes, PICRUSt could produce very wrong results.  You can use ids other than greengenes, but in that case you will be doing a custom analysis and will not be able to rely on our precalculated set of gene content predictions.  

* **Sequences that fail to map to references will not be predicted**: It is useful to check how well reference-based OTU picking worked.  Some samples, including for example many from marine environments, may not have representatives on the greengenes tree at a 97% OTU threshold.  So it is important to check what percentage of your reads were successfully mapped to greengenes.    If you use QIIME for your reference-based OTU picking step, this is available in the log file generated along with your OTU table.  If the most reads failed to map, you may wish to consider using a more generous similarity threshold (e.g. 94%).   It is important to be aware that this will have an effect on the accuracy of PICRUSt, since dissimilar organisms will be lumped together.  However, losing the vast majority of sequences may have an even greater effect.    The section on using the Nearest Sequenced Taxon Index to estimate overall PICRUSt accuracy will discuss how you can use that metric to estimate error rates.

3. Calculate reference genome coverage for your samples using NSTI scores
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When characterizing PICRUSt’s accuracy, we found that the largest single factor contributing to metagenome prediction accuracy for a given sample was the extent to which organisms from that sample had their genomes sequenced.  We developed the weighted Nearest Sequenced Taxon Index (weighted NSTI) score to summarize the extent to which microorganisms in a given sample are related to sequenced genomes.  NSTI scores are simply the average branch length that separates each OTU in your sample from a reference bacterial genome, weighted by the abundance of that OTU in the sample.   So the unit for the NSTI score is the same as your reference tree (typically 16S rRNA substitutions/site).   Thus a NSTI score of 0.03 indicates that the average microbe in your sample can be predicted using a relative from the same (97%) species.  

When comparing NSTI scores to typical bacterial taxonomic cutoffs (e.g. 97% sequence identity), it is worth noting that substitutions/site is equivalent to 1.0 - % identity at *short* phylogenetic distances, but greater than that value at longer phylogenetic distances due to, for example, multiple mutations at the same site that return the nucleotide to a previous state (Supp. Fig. 2 in `Zaneveld et al 2010 <http://nar.oxfordjournals.org/content/suppl/2010/03/02/gkq066.DC1/nar-02618-s-2009-File002.pdf>`_ has a figure relates these two for short distances) 

4. Calculate metagenomic confidence intervals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


A newer capability of PICRUSt is the ability to generate 95% confidence intervals for each gene prediction. These confidence intervals may be useful to examine for particular gene families of interest to determine the range of values that could be present in a community. 

**NOTE:** Although the code for prediction of metagenomic confidence intervals is unit-tested (produces correct output on small example datasets), and confidence intervals during generation of the precalculated gene content prediction files were benchmarked for accuracy against real genomes in the PICRUSt manuscript, metagenomic confidence intervals have not yet been biologically validated on large datasets.  Therefore this capability should be treated as experminetal for now, but may provide a rough guide for the range of values that might be present for each gene family in the metagenome.  See the :ref:`algorithm_description` for further details.

Metagenomic confidence intervals can be output using the --with_confidence option in the  predict_metagenomes.py script.


.. _empirical_qc_checklist:

Empirical Quality Control:  Sequencing
--------------------------------------

One direct method for testing PICRUSt metagenome prediction accuracy in a new environment is to sequence 16S rRNA and shotgun metagenomes from the same subset of  samples, and then directly compare the PICRUSt-predicted results against the empirical, sequenced results.  Of course this approach is only practical in cases where substantial collections of samples are slated for 16S rRNA sequencing.  But it can provide a great deal of empirical reassurance when sequencing many 16S samples, at the cost of perhaps a dozen metagenomes (which may be useful anyway for other purposes).  This approach was used to assess accuracy for the manuscript, and several PICRUSt scripts help to automate it.  

The most important script to check out is :ref:`compare_biom.py <compare_biom>`   This script allows you to compare how accurately one or more observed BIOM-format tables predict an expected table.  So you can run compare_biom.py using the BIOM table generated from your shotgun metagenome as the expected value, and your PICRUSt-predicted table as your observed result.   


Considerations:
     1.  **Control metagenomes should be fairly deep:**:  Because taxa saturate more quickly than genes, substantial metagenomic sequencing depth is needed to compare against PICRUSt predictions.  Based on rarefaction analysis of paired 16S rRNA/metagenome libraries from diverse soils, we find in the paper that roughly 72,000 raw or 15,000 annotated metagenomic sequences were needed before a subset of a deeply sequenced metagenome did better against the full metagenome than PICRUSt did (at least in soils).  That is, PICRUSt actually did better at predicting a deeply sequenced metagenome than did shallow subsets of the same metagenome. Therefore shallow metagenomes (although useful for assessing general functional changes in a community) make poor ‘gold standard’ positive controls for assessing PICRUSt accuracy.
     2.  **Randomized controls are needed:**: Even predicting random bacterial genomes produces substantial correlations-  simply knowing that a metagenome has mostly bacteria in it goes a long way towards correctly predicting the full metagenome.  This is simply because there are many common patterns in which gene families are more or less abundant (on average) across bacteria and archaea.  Therefore it is very useful to compare the accuracy of your PICRUSt prediction against the accuracy found when just predicting a random set of bacterial/archaeal genomes.   You can scramble the abundances of OTUs in your sample using the --random option in compare_biom.py.

.. _reporting_picrust_results:

Suggestions for Reporting Results From PICRUSt
----------------------------------------------

Inputs to PICRUSt -- picking reference OTUs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In addition to standard methods reporting on the methods used to sample, extract, amplify, and sequence your 16S rRNA data, some details of the reference-based OTU-picking step may be useful to report in order to ensure reproducibility of the results.  These include: 

    **The version of the Greengenes reference was used**:  You can run predict_metagenomes.py with the --help option to see a list of default parameters.  The Greengenes version is specified using the --gg_version option. If you didn't specifically request a particular version of Greengenes, the listed default value will be the version of greengenes used.
        
    **The proportion of reads that mapped to reference during OTU picking**: For some poorly explored environments, many or all sequences may fail to map to reference sequences.  Therefore it is a good idea to check the log files from QIIME's otu picking step to ensure that a reasonable number of sequences survived OTU picking.  Reporting these percentages may provide valuable information for other researchers trying to duplicate the analysis.

    **The similarity threshold used during OTU-picking**: typically 97% OTUs are used, 


PICRUSt-specific information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* What version of PICRUSt was used?  You can access this information by running:

  ``print_picrust_config.py``

* Were predictions corrected for predicted 16S rRNA copy number?  If you ran normalize_by_copy_number.py, then your predictions were normalized, and it will be useful to report this.

  * What was the Nearest Sequenced Taxon Index (NSTI) for the samples?   NSTI scores reflect the availability of reference genomes that are closely related to the most abundant microorganisms in your sample.  High scores (~ >0.15) generally mean few related references are available and predictions will be of low quality.  Low scores (roughly <0.06) indicate availability of closely related reference genomes.  NSTI scores can be calculated during metagenome prediction by passing the --with_accuracy option (see the :ref:`predict_metagenomes.py <predict_metagenomes>` script) 

* If you discuss a particular gene/KO/COG, you may be interested in the confidence interval for that metagenome prediction.  You can generate this information during the metagenome prediction step (see the :ref:`predict_metagenomes.py <predict_metagenomes>` script).





