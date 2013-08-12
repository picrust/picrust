.. include:: global.rst
.. _install:

Installing PICRUSt
==================

.. note:: Most users will not need to install PICRUSt, but can instead use the `online Galaxy version <http://huttenhower.sph.harvard.edu/galaxy/root?tool_id=PICRUSt_normalize>`_.

PICRUSt is written in python, and has been tested on Mac OS X and Linux systems. To install PICRUSt, first install all of the mandatory requirements, following the instructions on their respective websites. You should then download PICRUSt. You have the choice of downloading either the release version (recommended for most users) or the development version (recommended for users who need access to the latest features, and are willing to tolerate some instability in the code). Next, you will download the large precalculated PICRUSt files and place them in your picrust/data directory. Finally, you'll install PICRUSt. Each of these steps are detailed below.

Step 1. Install Requirements
----------------------------

Follow the install instructions found on the website of each of the dependencies below to install PICRUSt's dependencies.

**Mandatory**

* `python`_ (version 2.7)
* `PyCogent`_ (version 1.5.3)
* `biom`_ (version 1.1.2)
* `numpy`_ (version 1.5.1)

**Rebuilding PICRUSt OR Genome Prediction (optional)**

* `R`_ installed with `APE`_ library


Step 2. Download PICRUSt
------------------------

Release software
^^^^^^^^^^^^^^^^

The latest release of PICRUSt is 0.9.2. You can download this file from `this link <https://github.com/picrust/picrust/releases/download/0.9.2/picrust-0.9.2.tar.gz>`_. 

We recommend the release version of PICRUSt for most users. If you're not sure whether you want the release or the development version of PICRUSt, you should likely go with the release version.

Development software
^^^^^^^^^^^^^^^^^^^^

Alternatively you can download the latest development version of PICRUSt. You can download the development version of PICRUSt from `this link <https://github.com/picrust/picrust/archive/master.zip>`_, or using the following command if you have ``git`` installed::

	git clone git://github.com/picrust/picrust.git picrust

Step 3. Download PICRUSt's precalculated files
----------------------------------------------

PICRUSt precomputes most of the computationally intensive pipeline so each user can get predictions with less steps. These files are changed whenever there is a new release of the GreenGenes tree and/or a new release of IMG. These files are fairly large and need to be downloaded separately. Also, you only need the version that corresponds to the GreenGenes that you picked OTUs against (see :ref:`otu_picking_tutorial`).

Assuming you will be picking OTUs against the newest version of GreenGenes and want KEGG Ortholog predictions, then the following files will be sufficient. If you picked against an older version of GreenGenes or are interested in different functional predictions (e.g. COGs) then see the complete list of :ref:`picrust_precalculated_files`.

**These files must be placed in your `picrust-0.9.2/picrust/data` directory before installing.**

* Download for `16S copy number normalization <https://github.com/picrust/picrust/releases/download/0.9.2/16S_13_5_precalculated.tab.gz>`_

* Download for `KEGG Ortholog predictions <https://dl.dropbox.com/s/c4yntzmdbnq7f9n/ko_13_5_precalculated.tab.gz>`_ 


Step 4. Install PICRUSt
-----------------------

After downloading PICRUSt, you'll need to unzip the file. If you've downloaded the release version, do this with the following command::
	
	tar -xzf picrust-0.9.2.tar.gz

You'll then change into the new ``picrust-0.9.2`` directory as follows::
	
	cd picrust-0.9.2

And finally, you'll install PICRUSt with the following command::
	
	sudo python setup.py install

If you don't have ``sudo`` access on the system where you're trying to install PICRUSt, there are several variations on this command that you can use to install only for the current user (as opposed to a system-wide installation, which the above command will perform). See `this discussion <http://docs.python.org/2/install/index.html#alternate-installation>`_.

