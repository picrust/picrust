#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2012, The PICRUST project"
__credits__ = ["Greg Caporaso","Morgan Langille"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from numpy import dot, array
from biom.table import table_factory

def predict_metagenomes(otu_table,genome_table):
    """ predict metagenomes from otu table and genome table 
    """
    # identify the overlapping otus that can be used to predict metagenomes 
    overlapping_otus = list(set(otu_table.ObservationIds) & 
                            set(genome_table.ObservationIds))
    
    if len(overlapping_otus) < 1:
        raise ValueError,\
         "No common OTUs between the otu table and the genome table, so can't predict metagenome."
    
    # create lists to contain filtered data - we're going to need the data in 
    # numpy arrays, so it makes sense to compute this way rather than filtering
    # the tables
    otu_data = []
    genome_data = []
    
    # build lists of filtered data
    for obs_id in overlapping_otus:
        otu_data.append(otu_table.observationData(obs_id))
        genome_data.append(genome_table.observationData(obs_id))
    
    # matrix multiplication to get the predicted metagenomes
    new_data = dot(array(otu_data).T,array(genome_data)).T
    
    # return the result as a sparse biom table - the sample ids are now the 
    # sample ids from the otu table, and the observation ids are now the 
    # functions (i.e., samples) from the genome table
    return table_factory(new_data,otu_table.SampleIds,genome_table.SampleIds)
