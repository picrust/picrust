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

* `16S <ftp://ftp.microbio.me/pub/picrust-references/picrust-1.0.0/16S_13_5_precalculated.tab.gz>`__
* `KO <ftp://ftp.microbio.me/pub/picrust-references/picrust-1.0.0/ko_13_5_precalculated.tab.gz>`__
* `COG <ftp://ftp.microbio.me/pub/picrust-references/picrust-1.0.0/cog_13_5_precalculated.tab.gz>`__
* `RFAM <ftp://ftp.microbio.me/pub/picrust-references/picrust-1.0.0/rfam_13_5_precalculated.tab.gz>`__

If you would also like to get confidence intervals you will need these additional files:

* `KO CIs <ftp://ftp.microbio.me/pub/picrust-references/picrust-1.0.0/ko_13_5_precalculated_variances.tab.gz>`__
* `COG CIs <ftp://ftp.microbio.me/pub/picrust-references/picrust-1.0.0/cog_13_5_precalculated_variances.tab.gz>`__
* `RFAM CIs <ftp://ftp.microbio.me/pub/picrust-references/picrust-1.0.0/rfam_13_5_precalculated_variances.tab.gz>`__


Greengenes 18may2012
--------------------

* `16S <https://github.com/picrust/picrust/releases/download/0.9.2/16S_18may2012_precalculated.tab.gz>`__
* `KO <https://github.com/picrust/picrust/releases/download/0.9.2/ko_18may2012_precalculated.tab.gz>`__
* `COG <https://github.com/picrust/picrust/releases/download/0.9.2/cog_18may2012_precalculated.tab.gz>`__


