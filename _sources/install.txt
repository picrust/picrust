.. include:: global.rst
.. _install:

Installing PICRUSt
==================

.. note:: Most users will not need to install PICRUSt, but can instead use the `online Galaxy version <http://huttenhower.sph.harvard.edu/galaxy/root?tool_id=PICRUSt_normalize>`_.

PICRUSt is written in python, and has been tested on Mac OS X and Linux systems. To install PICRUSt, first install all of the mandatory requirements, following the instructions on their respective websites. You should then download PICRUSt. You have the choice of downloading either the release version (recommended for most users) or the development version (recommended for users who need access to the latest features, and are willing to tolerate some instability in the code). Finally, you'll install PICRUSt. Each of these steps are detailed below.

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

The latest release of PICRUSt is 0.9.1. You can download this file from :download:`this link <releases/picrust-0.9.1.tar.gz>`. 

We recommend the release version of PICRUSt for most users. If you're not sure whether you want the release or the development version of PICRUSt, you should likely go with the release version.

Development software
^^^^^^^^^^^^^^^^^^^^

Alternatively you can download the latest development version of PICRUSt. You can download the development version of PICRUSt from `this link <https://github.com/picrust/picrust/archive/master.zip>`_, or using the following command if you have ``git`` installed::

	git clone git://github.com/picrust/picrust.git picrust

Step 3. Install PICRUSt
-----------------------

After downloading PICRUSt, you'll need to unzip the file. If you've downloaded the release version, do this with the following command::
	
	tar -xzf picrust-0.9.1.tar.gz

You'll then change into the new ``picrust-0.9.1`` directory as follows::
	
	cd picrust-0.9.1

And finally, you'll install PICRUSt with the following command::
	
	sudo python setup.py install

If you don't have ``sudo`` access on the system where you're trying to install PICRUSt, there are several variations on this command that you can use to install only for the current user (as opposed to a system-wide installation, which the above command will perform). See `this discussion <http://docs.python.org/2/install/index.html#alternate-installation>`_.

