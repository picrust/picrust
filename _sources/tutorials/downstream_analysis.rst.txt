.. _downstream_analysis_guide:

Analyzing PICRUSt predicted metagenomes
=======================================

.. include:: ../global.rst

Once you have your PICRUSt metagenome predictions (see :ref:`metagenome_prediction_tutorial`), you can analyze your predicted metagenomes as you would analyze actual metagenomes. 

Here are some recommendations using both PICRUSt commands and other external tools:

Collapse predictions into pathways
----------------------------------
* Collapse the thousands of predicted functions into higher categories (e.g KOs into KEGG Pathways).
* See PICRUSt's :ref:`categorize_by_function.py <categorize_by_function>`.

Determine which OTUs are contributing to particular functions
-------------------------------------------------------------
* Identify which OTUs are contributing which functions .
* See PICRUSt's :ref:`metagenome_contributions.py <metagenome_contributions>`.

Analyze with QIIME
------------------
* Compute alpha diversity, beta diversity, generate various plots, and much more using `QIIME <http://www.qiime.org>`_.
* See :ref:`qiime_tutorial`.

Analyze with STAMP
------------------
* Browse results, create PCA and bar plots, and make statistical inferences between pairs of samples or multiple groups all within a graphical interface using `STAMP <http://kiwi.cs.dal.ca/Software/STAMP>`_.
* See :ref:`stamp_tutorial`. 

Analyzing metagenomes with HUManN, LEfSe, and GraPhlAn
------------------------------------------------------
* See :ref:`HUManN_tutorial`
