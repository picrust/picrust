.. _downstream_analysis_guide:

Analyzing PICRUSt predicted metagenomes
=======================================

.. include:: ../global.rst

Once you have your PICRUSt metagenome predictions (see :ref:`metagenome_prediction_tutorial`), you can analyze your predicted metagenomes as you would analyze actual metagenomes. 

Here are some recommendations using both PICRUSt commands and other external tools:

1) Collapse the thousands of predicted functions into higher categories using PICRUSt's :ref:`categorize_by_function.py <categorize_by_function>`.

2) Identify which OTUs are contributing which functions using PICRUSt's :ref:`metagenome_contributions.py <metagenome_contributions>`.

3) Analyze the metagenome predictions using `QIIME <http://www.qiime.org>`_ (see :ref:`qiime_tutorial`).

4) Browse results, create PCA and bar plots, and make statistical inferences between pairs of samples or multiple groups (defined by your sample metadata) using `STAMP <http://kiwi.cs.dal.ca/Software/STAMP>`_ (see :ref:`stamp_tutorial`). 


