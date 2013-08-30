.. _quality_control:

Quality Control of PICRUST Predictions
======================================

.. include:: ../global.rst

This section covers steps that can be taken to ensure PICRUSt predictions will be as accurate as possible, and to characterize how well or poorly the gene content of a given set of 16S rRNA samples can be predicted.

Overview
--------

PICRUSt’s main use is in estimating the bacterial and archaeal genes present in a microbial community metagenome, using 16S rRNA data.   The output is a table of gene family abundances for each sample, with the exact scheme used to delineate gene families chosen by you (options include KEGG, COG, PFAM and RFAM).  The accuracy of PICRUSt’s predicted metagenomes can vary quite a bit depending on several factors including, most importantly, the extent to which sequenced genomes are available for the most abundant species in your community of interest.   

More detail about the types of tests we have done on the algorithm, and the detailed results are available in the `PICRUSt manuscript in Nature Biotechnology <http://www.nature.com/nbt/journal/vaop/ncurrent/abs/nbt.2676.html>`_

Here, we will cover the most important conclusions from that research, and delve more into some extra steps you can take in order to predict how accurate PICRUSt will be for your samples of interest.   We have arranged this information as checklists:

The ‘basic’ checklists include tests that you can run along with your PICRUSt analysis without doing any additional sequencing.  

The ‘advanced’ checklist is intended for larger projects (with for example hundreds or thousands of 16S rRNA libraries and perhaps a dozen or more shotgun metagenomes available), and includes validations similar to those used in the PICRUSt manuscript.   At the bottom of the page, we include a checklist of items you may wish to include in your methods section when describing PICRUSt so that others can duplicate your findings.

Both the ‘basic’ and ‘advanced’ checklists assume you are doing the main PICRUSt workflow- predicting a ‘virtual metagenome’ from 16S rRNA data.  However, PICRUSt can in theory be used to predict any trait (continuous evolutionary character) across organisms.   Additional discussion of quality control steps you might take if you are using PICRUSt in a customized workflow are discussed in the ‘Custom trait prediction with PICRUSt’ section.

We would like this page to be a resource for the community.  If you have additional suggestions for additional steps that will help the community with their PICRUSt analysis, we would love to hear them-  please e-mail picrust-users@googlegroups.com

General limitations of PICRUSt
------------------------------

* Only genes from organisms hit by your primer will be included.  The main use for PICRUSt is in taking 16S rRNA data, and predicting a metagenome from that data using evolutionary modelling of how gene content has changed relative to sequenced genomes.   Therefore PICRUSt can only predict the portion of the metagenome that is contributed by the set of organisms picked up by your primers.   The PICRUSt validation datasets used universal 515f/806r V4 16S rRNA primers designed to minimize (though not eliminate) bias across bacterial/archaeal taxonomy [link Caporaso et al 2011, Walters et al 2011].   If the primers used don’t amplify an organism, then of course that organisms’ contribution to the metagenome is not predicted.   As an example, many popular 16S rRNA primers were did not work well for amplifying *Verrucomicrobia*.   If your sample had a large proporition of *Verrucomicrobia* in it, and you used such a primer set, then of course the metagenome predicted by PICRUSt would also fail to include genes contributed by *Verrucomicrobia*.   Similarly, since the input data for the standard workflow is 16S rRNA, any eukaryotic or viral contributions to the metagenome will not be predicted.

* PICRUSt can only predict gene families that are part of the input data.  Therefore, genes that do not fall within the orthology scheme used will not be predicted.  

* PICRUSt functional output is purely based on the particular input reference used. Therefore, any gaps or inaccuracies in pathway annotation or assignments of gene function will still be present.  As an example, many KEGG Orthology groups are listed as participating in pathways not found in bacteria or otherwise not reflective of true function.  In many cases this is simply due to bacteria containing  (distant) homologs of enzymes with important roles in mammalian pathways.  Therefore, it is worth carefully checking KEGG pathway annotations to ensure that they are reasonable for your system. 


Basic Quality Control steps for PICRUSt
---------------------------------------

These are some steps that can be taken to ensure that you get the most out of your PICRUSt analysis. 

* *Description of the Algorithm and Test Data* The PICRUSt algorithm and a number of accuracy tests are described in the :ref:`paper <http://www.nature.com/nbt/journal/vaop/ncurrent/abs/nbt.2676.html>`_. An additional description of the algorithm is available on the PICRUSt webpage here: :ref:`algorithm_description`.  Importantly, testing indicates that while PICRUSt can work quite well when adequete reference genomes are available, results may be unreliable for communities where adequete genome cover but also others where it is not reliable. 

Make sure you’ve done appropriate QC on the input 16S rRNA data.   Unless you are doing a customized analysis [link], PICRUSt requires that input 16S rRNA data be in the form of a BIOM [link] table that was picked by mapping your reads to references in the greengenes tree.   The most straightforward way to do this is using QIIME’s reference-based OTU picking pipeline as described here [link].   

 This has several implications:
        --OTU IDs must be based on greengenes: If you input an OTU table generated some other way, but with some IDs that happen to match greengenes, PICRUSt could produce very wrong results.  You can use ids other than greengenes, but in that case you will be doing a custom analysis[link] and will not be able to rely on our precalculated set of gene content predictions.  
        -- It is useful to check how well reference-based OTU picking worked.  Some samples, including for example many from marine environments, may not have representatives on the greengenes tree at a 97% OTU threshold.  So it is important to check what percentage of your reads were successfully mapped to greengenes.    If you use QIIME for your reference-based OTU picking step, this is available in the log file generated along with your OTU table.  If the most reads failed to map, you may wish to consider using a more generous similarity threshold (e.g. 94%).   It is important to be aware that this will have an effect on the accuracy of PICRUSt, since dissimilar organisms will be lumped together.  However, losing the vast majority of sequences may have an even greater effect.    The section on using the Nearest Sequenced Taxon Index to estimate overall PICRUSt accuracy will discuss how you can use that metric to estimate error rates.
       -- Calculate NSTI scores to check how closely the dominant microorganisms in your samples are related to sequenced bacterial or archaeal genomes.     When characterizing PICRUSt’s accuracy, we found that the largest single factor contributing to metagenome prediction accuracy for a given sample was the extent to which organisms from that sample had their genomes sequenced.  We developed the weighted Nearest Sequenced Taxon Index (weighted NSTI) score to summarize the extent to which microorganisms in a given sample are related to sequenced genomes.  NSTI scores are simply the average branch length that separates each OTU in your sample from a reference bacterial genome, weighted by the abundance of that OTU in the sample.   So the unit for the NSTI score is the same as your reference tree, and is typically 16S rRNA substitutions per site.   Thus a NSTI score of 0.03 indicates that the average microbe in your sample can be predicted using a relative from the same (97%) species

NSTI Score
Spearman correlation: predicted vs. expected

Advanced Quality Control:  Sequencing
-------------------------------------

In our view, the most direct and reliable way to test whether PICRUSt can predict metagenomes in a new environment is to sequence 16S rRNA and shotgun metagenomes from the same samples, and then directly compare the PICRUSt-predicted results against the empirical, sequenced results.  Of course this approach is only possible if you plan on doing substantial 16S rRNA sequencing.  But it can provide a great deal of empirical reassurance when sequencing many 16S samples, at the cost of perhaps a dozen metagenomes (which may be useful anyway for other purposes).  This approach was used to assess accuracy for the manuscript, and several picrust scripts help to automate it.  

The most important script to check out is compare_biom.py   This script allows you to compare how accurately one or more observed BIOM-format tables predict an expected table.  So you can run compare_biom.py using the BIOM table generated from your shotgun metagenome as the expected value, and your PICRUSt-predicted table as your observed result.   


Considerations:
     1.  Control metagenomes should be fairly deep:  Because taxa saturate more quickly than genes, you need substantial metagenomic sequencing depth to compare against PICRUSt predictions.  Based on rarefaction analysis of paired 16S rRNA/metagenome libraries from diverse soils, we find in the paper that for that environment roughly 72,000 raw or 15,000 annotated metagenomic sequences were needed before a subset of a deeply sequenced metagenome did better against the full metagenome than PICRUSt did.  That is, PICRUSt actually did better at predicting a deeply sequenced metagenome than did shallow subsets of the same metagenome. Therefore shallow metagenomes (although useful for assessing general functional changes in a community) make poor ‘gold standard’ positive controls for assessing PICRUSt accuracy.
     2.  Even predicting random bacterial genomes produces substantial correlations:  Simply knowing that a metagenome has mostly bacteria in it goes a long way towards correctly predicting the full metagenome.  This is simply because there are many common patterns in which gene families are more or less abundant (on average) across bacteria and archaea.  Therefore it is very useful to compare the accuracy of your PICRUSt prediction against the accuracy found when just predicting a random set of bacterial/archaeal genomes.   You can scramble the abundances of OTUs in your sample using the --random option in compare_biom.py.

Suggestions for Reporting Results From PICRUSt
------------------------------

General
^^^^^^^
* Standard reporting on the methods used to sample, extract, amplify, and sequence your 16S rRNA data.
* The steps used in processing your sequenced reads, and mapping them to a reference.
    Per the discussion above, what proportion of reads mapped to reference?  What similarity threshold was used?

PICRUSt-specific information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* What version of PICRUSt and greengenes was used?  You can access this information by running print_picrust_config.py [check]
* What was the Nearest Sequenced Taxon Index (NSTI) for your samples?   NSTI scores reflect the availability of reference genomes that are closely related to the most abundant microorganisms in your sample.  High scores (~ >0.15) generally mean few related references are available and predictions will be of low quality.  Low scores (roughly <0.06) indicate availability of closely related reference genomes.  
* If you discuss a particular gene/KO/cog, you may be interested in the confidence interval for that metagenome prediction.  You can generate this information during the metagenome prediction step (see :ref:`predict_metagenomes`).
* Were predictions corrected for predicted 16S rRNA copy number?  If you ran normalize_by_copy_number.py, then your predictions were normalized, and it will be useful to report this




