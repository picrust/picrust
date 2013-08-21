.. _stamp_tutorial:

Analyzing metagenomes with STAMP
================================

.. include:: ../global.rst

`STAMP <http://kiwi.cs.dal.ca/Software/STAMP>`_ (Statistical Analysis of Metagenomic Profiles) is a software package for analyzing metagenomic profiles. A user friendly, graphical interface permits easy exploration of statistical results and generation of publication quality plots for inferring the biological relevance of features in a metagenomic profile.

The metagenome predictions from PICRUSt, can be used with slight altering of the tab-delimited output files (using the ``-f`` option in predict_metagenomes.py and categorize_by_function.py). 

STAMP takes as input two data files. The first is a ``profile file`` that contains the count data, and the second is an optional ``group metadata file`` that contains metadata about each of the samples to allow comparison of groups (e.g. sick vs healthy). 

The ``group metadata file`` will need to be created by yourself, but is a simple tab-delimited format that can be easily exported from a spreadsheet program. The first column must have sample ids that match the ``profile file`` with each of the other columns representing different metadata about the samples.

An example ``group metadata file`` could look like this::

    Sample Id    Location    Phenotype    Gender    Sample Size
    Sample 1     Canada        Obese      Female      4000
    Sample 2     Canada        Lean       Male        2000
    Sample 3     Italy         Lean       Female      3000 

To create the ``profile file`` you need to do a simple conversion of PICRUSt's output. The first header line should be discarded and the last column containing metadata should be discarded. This can be done in a simple text editor and spreasheet program or you can use the unix command below.

Assuming your output file was created with the following command::

    predict_metagenomes.py -f -i normalized_otus.biom -o predicted_metagenomes.txt

then this file can be converted to a STAMP ``profile file``::

    sed '1d' predicted_metagenomes.txt | rev | cut -f 2- | rev > predicted_metagenome.spf

This will allow STAMP analysis on each of the ~7000 KOs. However, it may be beneficial to first collapse the counts into more general categories::

    categorize_by_function.py -f -i predicted_metagenomes.biom -c KEGG_Pathways -l 3 -o predicted_metagenomes.L3.txt

then convert this file to STAMP ``profile file``::

    sed '1d' predicted_metagenomes.L3.txt | rev | cut -f 2- | rev > predicted_metagenome.L3.spf

Once you have your STAMP ``profile file`` and ``group metadata file`` generated, then install STAMP, load the files from the **File->Load Data...** menu and start browsing your data. 

The `STAMP User's Guide <http://kiwi.cs.dal.ca/wow/images/0/02/STAMP_Users_Guide_v2.0.0.pdf>`_ is very clear and helpful. 

With only a few clicks you can generate figures like these:

.. image:: /images/STAMP2_Windows7_Screenshot.png
	   :width: 600px
	   
.. image:: /images/ArumugamBoxPlot.png
	   :width: 600px
	   
.. image:: /images/STAMP2_RumenErrorBarPlot.png
	   :width: 600px
	   

