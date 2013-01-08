#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Greg Caporaso","Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
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

             
def extract_otu_and_genome_data(otu_table,genome_table):
    """Return lists of otu,genome data, and overlapping genome/otu ids
    
    otu_table -- biom Table object for the OTUs
    genome_table -- biom Table object for the genomes
    """
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
    return otu_data,genome_data,overlapping_otus 

        

def predict_metagenomes(otu_table,genome_table,verbose=False):
    """ predict metagenomes from otu table and genome table 
    """
    
    otu_data,genome_data,overlapping_otus = extract_otu_and_genome_data(otu_table,genome_table)
    # matrix multiplication to get the predicted metagenomes
    new_data = dot(array(otu_data).T,array(genome_data)).T
    
    #Round counts to nearest whole numbers
    new_data = around(new_data)
    
    # return the result as a sparse biom table - the sample ids are now the 
    # sample ids from the otu table, and the observation ids are now the 
    # functions (i.e., observations) from the genome table


    result_table =  table_factory(new_data,otu_table.SampleIds,genome_table.ObservationIds)

    #We need to preserve metadata about the samples from the OTU table, 
    #and metadata about the gene functions from the genome table
    
    #Transfer sample metadata from the OTU table
    #to the metagenome table (samples are the same)
    result_table = transfer_metadata(otu_table,result_table,\
      donor_metadata_type='SampleMetadata',\
      recipient_metadata_type='SampleMetadata',verbose=verbose)
    
    #Now transfer observation metadata (e.g. gene metadata) 
    #from the genome table to the result table
    result_table = transfer_metadata(genome_table,result_table,\
      donor_metadata_type='ObservationMetadata',\
      recipient_metadata_type='ObservationMetadata',verbose=verbose)
    

    return result_table


def transfer_metadata(donor_table,recipient_table,\
        donor_metadata_type="ObservationMetadata",recipient_metadata_type="ObservationMetadata",\
        verbose = True):
    """Transfer particular metadata properties from donor table to recipient BIOM table"""
    #Check that property label looks OK
    for metadata_type in [donor_metadata_type,recipient_metadata_type]:
        if metadata_type not in ['ObservationMetadata','SampleMetadata']:
            raise ValueError("Donor and recipient metadata types passed to transfer_metadata must be either 'ObservationMetadata' or 'SampleMetadata'")
   
    if verbose: 
        print "Transferring donor_table.%s to recipient_table.%s" %(donor_metadata_type,recipient_metadata_type)
    
    if donor_metadata_type == "ObservationMetadata":
        transfer_fn = transfer_observation_metadata
    elif donor_metadata_type == "SampleMetadata":
        transfer_fn = transfer_sample_metadata
    
    recipient_table = transfer_fn(donor_table,recipient_table,recipient_metadata_type=recipient_metadata_type,verbose=verbose)
    return recipient_table

def transfer_observation_metadata(donor_table,recipient_table,\
        recipient_metadata_type="ObservationMetadata",\
        verbose = True):
    """Transfer observation metadata properties from donor table to recipient BIOM table"""
    #Check that property label looks OK
    donor_metadata_type = "ObservationMetadata"
    if recipient_metadata_type not in ['ObservationMetadata','SampleMetadata']:
            raise ValueError("Recipient metadata type passed to transfer_metadata must be either 'ObservationMetadata' or 'SampleMetadata'")
   
    if verbose: 
        print "Transferring donor_table.%s to recipient_table.%s" %(donor_metadata_type,recipient_metadata_type)
    
    donor_metadata = getattr(donor_table,donor_metadata_type,None) 
   
    if not donor_metadata:
        #No metadata to transfer, so nothing more needs to be done.
        return recipient_table

    metadata = {}
    md_ids = donor_table.ObservationIds

    for md_id in md_ids:
        metadata_value = donor_table.ObservationMetadata[donor_table.getObservationIndex(md_id)]
        metadata[str(md_id)] = metadata_value
    if recipient_metadata_type == "ObservationMetadata":
        recipient_table.addObservationMetadata(metadata)
    elif recipient_metadata_type == "SampleMetadata":
        recipient_table.addSampleMetadata(metadata)
   
    return recipient_table

def transfer_sample_metadata(donor_table,recipient_table,\
        recipient_metadata_type="SampleMetadata",\
        verbose = True):
    """Transfer sample metadata properties from donor table to recipient BIOM table"""
    #Check that property label looks OK
    donor_metadata_type = "SampleMetadata"
    if recipient_metadata_type not in ['ObservationMetadata','SampleMetadata']:
            raise ValueError("Recipient metadata type passed to transfer_metadata must be either 'ObservationMetadata' or 'SampleMetadata'")
   
    if verbose: 
        print "Transferring donor_table.%s to recipient_table.%s" %(donor_metadata_type,recipient_metadata_type)
    
    donor_metadata = getattr(donor_table,donor_metadata_type,None) 
   
    if not donor_metadata:
        #No metadata to transfer, so nothing more needs to be done.
        return recipient_table

    metadata = {}
    md_ids = donor_table.SampleIds

    for md_id in md_ids:
        metadata_value = donor_table.SampleMetadata[donor_table.getSampleIndex(md_id)]
        metadata[str(md_id)] = metadata_value
    
    if recipient_metadata_type == "SampleMetadata":
        recipient_table.addSampleMetadata(metadata)
    elif recipient_metadata_type == "ObservationMetadata":
        recipient_table.addObservationMetadata(metadata)
   
    return recipient_table


def calc_nsti(otu_table,genome_table,weighted=True):
    """Calculate the weighted Nearest Sequenced Taxon Index for ids
    otu_table -- a .biom OTU table object
    genome_table -- the corresponding set of PICRUST per-OTU genome predictions
    weighted -- if True, normalize by OTU abundance
    """
    
    # identify the overlapping otus that can be used to calculate the NSTI
    overlapping_otus = get_overlapping_ids(otu_table,genome_table)
    total = 0.0 
    n = 0.0
    observation_ids = otu_table.SampleIds
    for obs_id in overlapping_otus:
        #print "Curr observed id:",obs_id
        obs_id_idx = genome_table.getSampleIndex(obs_id)
        #observatin_ids.append(genome_table.ObservationIds[obs_id_idx])
        curr_nsti =  genome_table.SampleMetadata[obs_id_idx]['NSTI']
        #print "Current NSTI", curr_nsti
        if weighted:
            curr_counts = otu_table.observationData(obs_id)
            #print "Curr counts:",curr_counts
            total += curr_counts*curr_nsti
            #print "total:",total
            n += curr_counts
            #print "n:",n
        else:
            total += curr_nsti
            n += 1
    
    result=total/n
    #print "result:",result
    
    return observation_ids,result

