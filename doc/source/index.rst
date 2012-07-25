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

Download
========

Two files are needed to **run PICRUST** on your own data:

* The `PICRUST software`_ is available from our `Sourceforge page <http://picrust.sourceforge.net>`_
* To generate a PICRUST-compatible OTU table, you need this `reference data for OTU picking  <https://s3.amazonaws.com/picrust-public-data/img_gg_otus_18may2012.tgz>`_, the use of which is described in the :ref:`otu_picking_tutorial`.

Additional files are available for those looking to **develop for PICRUST** or perform **single-genome ancestral state reconstruction**:

* PICRUST uses a suite of `ancestral state reconstruction tools <http://***>`_ to generate the data needed for metagenomic inference.

Quickstart
==========

Looking to get metagenome predictions for your 16S data? Follow the :ref:`quickstart`.

Contact Us
==========

For PICRUST announcements and questions, including notification of new releases, you can subscribe to the `PICRUST users list <https://groups.google.com/group/picrust-users/subscribe?note=1&hl=en&noredirect=true&pli=1>`_.

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
