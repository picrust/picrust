.. include:: global.rst

PICRUST: Phylogenetic Investigation of Communities by Reconstruction of Unobserved STates
=========================================================================================

PICRUST is a bioinformatics software package, **currently still in development**, designed to predict metagenome functional content from marker gene (e.g., 16S rRNA) surveys and full genomes.

PICRUST is freely available here under the `GPL`_. However, the software and documentation are still under active development, so we discourage the use of PICRUST for general scientific research (for now).

.. toctree::
   :maxdepth: 2
   
   ./tutorials/index
   ./methods/index
   ./scripts/index

Requirements
============
**Mandatory**

* `PyCogent`_
* `biom`_

**Rebuilding PICRUST OR Genome Prediction (optional)**

* `R`_ installed with `APE`_ library


Download
========

**PICRUST Software (All Users)**

* PICRUST software

.. warning::

       An official release version is not yet available. Please instead download the `PICRUST development software`_.

**OTU table preperation (Most Users)**

* Before using PICRUST for metagenome prediction, you must ensure your OTU table is PICRUST-compatible (OTU identifiers are `Greengenes`_ identifiers). 
* To generate a PICRUST-compatible OTU table, you need the `PICRUST GG reference data`_. 
* For use, see tutorial: :ref:`otu_picking_tutorial`

**Inferring Individual Genomes (Optional: For Advanced Users)**

* Additional files are available for those looking to **develop for PICRUST** or perform **single-genome ancestral state reconstruction**: `PICRUST starting files`_.
* For use, see tutorial: :ref:`genome_prediction_tutorial`


Install
=======

* Add the location of your 'PICRUST' directory to your PYTHONPATH environment variable.
* Add the location of your 'PICRUST/scripts' to your PATH environment variable for ease of use (optional).

Quickstart
==========

Looking to get metagenome predictions for your 16S data? Follow the :ref:`quickstart_tutorial`.

Contact Us
==========

For PICRUST announcements and questions, including notification of new releases, you can subscribe to the `PICRUST users list`_.

Citation
========

Hopefully coming soon ;)

News & Announcements
====================

.. include:: news.rst
   :start-line: 4
   :end-line: 7
   
* More :ref:`news`

..	Indices and tables
	==================
	
	* :ref:`genindex`
	* :ref:`modindex`
	* :ref:`search`
