.. include:: global.rst
.. _picrust_precalculated_files:

PICRUSt's Precalculated Files
=============================

PICRUSt requires the downloading of precalculated files. 

These downloaded files should be placed within your `picrust/data` directory BEFORE installation.

Download the set of files that correspond to the version of Greengenes that you used for OTU picking (see :ref:`otu_picking_tutorial`)

The minimum set of files are:

1) 16S for copy number normalization (used in the script :ref:`normalize_by_copy_number.py <normalize_by_copy_number>`).
2) Whatever type of function predictions you want returned by PICRUSt (e.g. KO or COG). This is used in the script :ref:`predict_metagenomes.py <predict_metagenomes>`.


Greengenes v13.5 (and IMG 4)
----------------------------

* `16S <ftp://thebeast.colorado.edu/pub/picrust-references/picrust-0.9.2/16S_13_5_precalculated.tab.gz>`_
* `KO <ftp://thebeast.colorado.edu/pub/picrust-references/picrust-0.9.2/ko_13_5_precalculated.tab.gz>`_
* `COG <ftp://thebeast.colorado.edu/pub/picrust-references/picrust-0.9.2/cog_13_5_precalculated.tab.gz>`_
* `RFAM <ftp://thebeast.colorado.edu/pub/picrust-references/picrust-0.9.2/rfam_13_5_precalculated.tab.gz>`_

If you would also like to get confidence intervals you will need these additional files:

* `KO CIs <ftp://thebeast.colorado.edu/pub/picrust-references/picrust-0.9.2/ko_13_5_precalculated_variances.tab.gz>`_
* `COG CIs <ftp://thebeast.colorado.edu/pub/picrust-references/picrust-0.9.2/cog_13_5_precalculated_variances.tab.gz>`_
* `RFAM CIs <ftp://thebeast.colorado.edu/pub/picrust-references/picrust-0.9.2/rfam_13_5_precalculated_variances.tab.gz>`_


Greengenes 18may2012
--------------------

* `16S <https://github.com/picrust/picrust/releases/download/0.9.2/16S_18may2012_precalculated.tab.gz>`_
* `KO <https://github.com/picrust/picrust/releases/download/0.9.2/ko_18may2012_precalculated.tab.gz>`_
* `COG <https://github.com/picrust/picrust/releases/download/0.9.2/cog_18may2012_precalculated.tab.gz>`_


