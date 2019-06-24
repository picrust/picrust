.. include:: global.rst
.. _install:

Installing PICRUSt
==================

.. note:: Most users will not need to install PICRUSt, but can instead use the `online Galaxy version <http://galaxy.morganlangille.com/>`_.

PICRUSt is written in python, and has been tested on Mac OS X and Linux systems. The easist way to install PICRUSt is with `conda` as described below. You will need to download the pre-calculated files before running the tool no matter how it is installed.

Recommended: Install PICRUSt from bioconda
----------------------------
The easiest way to install PICRUSt is from `bioconda`_, which requires `conda` to be installed (see `miniconda`_).

To do so you need to type these commands to install PICRUSt in a new environment named `picrust1`::

        conda create -n picrust1 -c bioconda picrust
        conda activate picrust1

You will then need to download the required precalculated files with this command (see more details under Step 4 below)::

	download_picrust_files.py

You should be good to go! If you run into any problems you can try install from source manually instead.

Install From Source Step 1. Install Requirements
----------------------------

Follow the install instructions found on the website of each of the dependencies below to install PICRUSt's dependencies.

**Mandatory**

h5py needs to be installed manually (numpy is also required, but it should be installed automatically with this package):

* `h5py`_ (version 2.7.1)

You can install h5py with this command::
	
        pip install h5py 
	    	
**Rebuilding PICRUSt's precalculated 16S rRNA OR Genome Predictions (optional)**

* `R`_ installed with `APE`_ library


Install From Source Step 2. Download PICRUSt
------------------------

Release software
^^^^^^^^^^^^^^^^

The latest release of PICRUSt is `picrust-1.1.4 <https://github.com/picrust/picrust/releases/download/v1.1.4/picrust-1.1.4.tar.gz>`_. 

We recommend the release version of PICRUSt for most users. If you're not sure whether you want the release or the development version of PICRUSt, you should likely go with the release version.

Development software
^^^^^^^^^^^^^^^^^^^^

Alternatively you can download the latest development version of PICRUSt. You can download the development version of PICRUSt from `this link <https://github.com/picrust/picrust/archive/master.zip>`_, or using the following command if you have ``git`` installed::

	git clone git://github.com/picrust/picrust.git picrust

Install From Source Step 3. Install PICRUSt
-----------------------

After downloading PICRUSt, you'll need to unzip the file. If you've downloaded the release version, do this with the following command::

	tar -xzf picrust-1.1.4.tar.gz

You'll then change into the new ``picrust-1.1.4`` directory as follows::

	cd picrust-1.1.4

And finally, you'll install PICRUSt with the following command::

	pip install .

These dependencies are automatically installed with `pip install .`:

* `python`_ (version 2.7)
* `PyCogent`_ (version 1.5.3)
* `biom`_ (version 2.1.7)


Install from Source Step 4. Download PICRUSt's precalculated files
----------------------------------------------

PICRUSt precomputes most of the computationally intensive pipeline so each user can get predictions with less steps. These files are changed whenever there is a new release of the GreenGenes tree and/or a new release of IMG. These files are fairly large and need to be downloaded separately. Also, you only need the version that corresponds to the GreenGenes that you picked OTUs against (see :ref:`otu_picking_tutorial`).

Assuming you will be picking OTUs against the newest version of GreenGenes and want KEGG Ortholog predictions, then the default files will be sufficient. If you picked against an older version of GreenGenes or are interested in different functional predictions (e.g. COGs) then see the complete list of :ref:`picrust_precalculated_files`.

You can install the default files automatically by running::

	download_picrust_files.py

Alternatively, download these files and place them in the ``data`` directory where PICRUSt was installed (e.g. if installed in a conda enviroment named "picrust1", place the files in ``$HOME/miniconda3/envs/picrust1/lib/python2.7/site-packages/picrust/data/``):

* Download for `16S copy number normalization <http://kronos.pharmacology.dal.ca/public_files/picrust/picrust_precalculated_v1.1.4/13_5/16S_13_5_precalculated.tab.gz>`_

* Download for `KEGG Ortholog predictions <http://kronos.pharmacology.dal.ca/public_files/picrust/picrust_precalculated_v1.1.4/13_5/ko_13_5_precalculated.tab.gz>`_ 
