.. _josh:

=====================================
Pre-PI-CRUST from the Huttenhower Lab
=====================================

This document covers the code in the ``pre_picrust/huttenhower/`` directory. Questions on usage should be directed to Joshua Reyes (jreyes@post.harvard.edu).

Description of code
===================

This code performs the following functions that are the backbone of PI-CRUST:
 * Constructs a (NCBI) labeled phylogeny from the input OTU abundance table.
 * Performs ancestral state reconstruction from the tree using a decaying exponential weighted vote to produce a gene count for each missing genome.
 * Predicts functional abundances simply by multiplying OTU abundance by gene count and summing across all OTUs.

Software Requirements
=====================
 * Python >= 2.6
 * Numpy, a package for scientific computing in Python
 * ETE 2, a package for the manipulation, analysis, and visualization of phylogenetic trees
 * SCons, a software construction tool similar to the Make utility
 
Additionally, I like to use:

 * virtualenv, a tool to create isolated Python environments
 * virtualenvwrapper, syntactic sugar for virtualenv
 * pip, a package manager for Python
 
Installation
============
PICRUSt itself doesn't require an installation. Simply copy our code into a directory $PICRUST_HOME. To get things running, though, you need to install the packages PICRUSt depends on.

Step 1: Create an isolated Python environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Set up a Python environment specifically to run PICRUST. This keeps the packages you install for PICRUSt separate from all of your other projects::

   mkvirtualenv picrust


Step 2: Install the dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To install the dependencies, I like to use pip.

I've included a copy of SCons 2.1.0 and ETE2 2.0 rev 111 in the etc/ folder for your convenience. You can rely on pip to install Numpy.
::
  pip install numpy
  cd $PICRUST_HOME/etc
  tar xzfv scons-2.1.0.tar.gz
  tar xzfv ete2-2.0rev111.tar.gz
  cd scons-2.1.0
  python setup.py install
  cd ../ete2-2.0rev111
  python setup.py install
  
PICRUSt doesn't use the database or visualization ETE toolkits. The ETE setup script will warn you that you won't be able to use them unless you install more packages. That's okay, though. Just answer yes to continue with its installation.
::
  Installing ETE (A python Environment for Tree Exploration).

  Checking dependencies...
  MySQLdb cannot be found in your python installation.
  MySQLdb is required for the PhylomeDB access API.
  PyQt4 cannot be found in your python installation.
  PyQt4 is required for tree visualization and rendering.

  However, you can still install ETE without such funtionalites.
  Do you want to continue with the installation? [y,n] y

Step 3: Make ANFunc path accessible
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
PICRUSt relies on an a binary file called ANFunc. This executable isn't actually part of PIECRUST; it's part of another project in the lab. It was written in C++ to multiply matrices fast. It can probably be replaced with something in SciPy later on.

For example::

  echo "export PATH=$PICRUST_HOME/bin/:$PATH" >> ~/.bashrc

Running PICRUSt
===============
Input
^^^^^
Simply drop as many OTU abundance PCL files as you like into the $PICRUST/input/ directory. Two other files, a phylogeny of taxa ``$PICRUST/taxdump.txt`` and KO table ``genomes_02-zro-flt.txt`` should also be included. You can find the KO table in the Dropbox::

    cp utility_data_files/huttenhower_prepicrust_input/* pre_picrust/huttenhower/input

Running
^^^^^^^
SCons drives the entire pipeline, which makes for a simple execution. Once everything is in place in the input directory, run ``scons`` from $PICRUST_HOME. The ``SConstruct`` file directs ``scons`` and can be modified to suit your needs.
::
    scons

Output
^^^^^^
Final results are saved in ``output/<input name>_kos.pcl`` by default.

Intermediate calculations are saved in ``output`` folder. Trees are saved using an internal representation in ``.fir`` files. Intermediate tables have end in the suffix ``_fir.pcl``.
