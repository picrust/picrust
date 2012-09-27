#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2012, The PICRUST project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from numpy import dot, array, around
from biom.table import table_factory

def get_overlapping_ids(otu_table,genome_table):
    """Get the ids that overlap between the OTU and genome tables"""
    overlapping_otus = list(set(otu_table.ObservationIds) & 
                            set(genome_table.SampleIds))
    
    if len(overlapping_otus) < 1:
        raise ValueError,\
         "No common OTUs between the otu table and the genome table, so can't predict metagenome."
    return overlapping_otus

def predict_metagenomes(otu_table,genome_table):
    """ predict metagenomes from otu table and genome table 
    """
    # identify the overlapping otus that can be used to predict metagenomes 
    overlapping_otus = get_overlapping_ids(otu_table,genome_table)
    # create lists to contain filtered data - we're going to need the data in 
    # numpy arrays, so it makes sense to compute this way rather than filtering
    # the tables
    otu_data = []
    genome_data = []
    
    # build lists of filtered data
    for obs_id in overlapping_otus:
        otu_data.append(otu_table.observationData(obs_id))
        genome_data.append(genome_table.sampleData(obs_id))
    
    # matrix multiplication to get the predicted metagenomes
    new_data = dot(array(otu_data).T,array(genome_data)).T
    
    #Round counts to nearest whole numbers
    new_data = around(new_data)
    
    # return the result as a sparse biom table - the sample ids are now the 
    # sample ids from the otu table, and the observation ids are now the 
    # functions (i.e., observations) from the genome table
    return table_factory(new_data,otu_table.SampleIds,genome_table.ObservationIds)

def calc_nsti(otu_table,genome_table,weighted=True):
    """Calculate the weighted Nearest Sequenced Taxon Index for ids
    otu_table -- a .biom OTU table object
    genome_table -- the corresponding set of PICRUST per-OTU genome predictions
    weighted -- if True, normalize by OTU abundance
    """
    
    # identify the overlapping otus that can be used to calculate the NSTI
    overlapping_otus = get_overlapping_otus(otu_table,genome_table)
    total = 0.0 
    n = 0.0
    for obs_id in overlapping_otus:
        curr_nsti =  otu_table.ObservationMetadata[obs_id]['NSTI']
        if weighted:
            curr_counts = otu_table.observationData(obs_id)
            total += counts*curr_nsti
            n += counts
        else:
            total += curr_nsti
            n += 1
    
    result=total/n
    return result

