#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Greg Caporaso","Jesse Zaneveld","Morgan Langille"]
__license__ = "GPL"
__version__ = "0.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from numpy import dot, array, around, asarray, sum as numpy_sum,sqrt
from biom.table import table_factory,SparseGeneTable, DenseGeneTable
from biom.parse import parse_biom_table, get_axis_indices, direct_slice_data, direct_parse_key
from picrust.predict_traits import variance_of_weighted_mean,calc_confidence_interval_95

def get_overlapping_ids(otu_table,genome_table,genome_table_ids="SampleIds",\
  otu_table_ids="ObservationIds"):
    """Get the ids that overlap between the OTU and genome tables
    otu_table - a BIOM Table object representing the OTU table
    genome_table - a BIOM Table object for the genomes
    genome_table_ids - specifies whether the ids of interest are SampleIds or ObservationIds
      in the genome table.  Useful because metagenome tables are often represented with genes
      as observations.
    """
    for id_type in [genome_table_ids,otu_table_ids]:
        if id_type not in ["SampleIds","ObservationIds"]:
            raise ValueError(\
              "%s is not a valid id type. Choices are 'SampleIds' or 'ObservationIds'"%id_type)
    overlapping_otus = list(set(getattr(otu_table,otu_table_ids)) & 
                            set(getattr(genome_table,genome_table_ids)))
    
    if len(overlapping_otus) < 1:
        print "OTU Ids:",getattr(otu_table,otu_table_ids)
        print "Genome Ids:",getattr(genome_table,genome_table_ids)
        raise ValueError,\
         "No common OTUs between the otu table and the genome table, so can't predict metagenome."
    return overlapping_otus

             
def extract_otu_and_genome_data(otu_table,genome_table,genome_table_ids="SampleIds",\
  otu_table_ids="ObservationIds"):
    """Return lists of otu,genome data, and overlapping genome/otu ids
    
    otu_table -- biom Table object for the OTUs
    genome_table -- biom Table object for the genomes
    """
    overlapping_otus = get_overlapping_ids(otu_table,genome_table,\
      genome_table_ids=genome_table_ids,otu_table_ids=otu_table_ids)
    # create lists to contain filtered data - we're going to need the data in 
    # numpy arrays, so it makes sense to compute this way rather than filtering
    # the tables
    otu_data = []
    genome_data = []
    
    # build lists of filtered data
    for obs_id in overlapping_otus:
        if otu_data == "SampleIds":
           otu_data.append(otu_table.sampleData(obs_id))
        else:
            otu_data.append(otu_table.observationData(obs_id))
        if genome_table_ids == "ObservationIds":
           genome_data.append(genome_table.observationData(obs_id))
        else:
           genome_data.append(genome_table.sampleData(obs_id))
    return otu_data,genome_data,overlapping_otus 


def load_subset_from_biom_str(biom_str,ids_to_load,axis="samples"):
    """Load a biom table containing subset of samples or observations from a BIOM format JSON string"""
    if axis not in ['samples','observations']:
        raise InputError(\
          'load_subset_from_biom_str axis parameter must be either "samples" or "observations"')

    ids = map(str,[l.strip() for l in ids_to_load])

    idxs, new_axis_md = get_axis_indices(biom_str, ids, axis)
    new_data = direct_slice_data(biom_str, idxs, axis)


    new_table_pieces = yield_subset_biom_str(biom_str,new_data,new_axis_md,axis)
    subset_biom_str = ''.join(new_table_pieces)
    return parse_biom_table(subset_biom_str)

def yield_subset_biom_str(biom_str,new_data,new_axis_md,axis):
    """Sequentially yield the components of a reduced biom string"""

    if axis not in ['samples','observations']:
        raise InputError(\
          'yield_subset_biom_str axis parameter must be either "samples" or "observations"')

    # Still has multiple walks over the file. bad form, but easy right now
    yield "{"
    keys = ["id","format","format_url","type",\
      "generated_by","date","matrix_type",\
      "matrix_element_type"]

    for key in keys:
        yield direct_parse_key(biom_str,key)
        yield ","
    
    for entry in [new_data,new_axis_md]:
        yield entry
        yield ","

    if axis == "observations":
        yield direct_parse_key(biom_str, "columns")
    else:
        yield direct_parse_key(biom_str, "rows")
    yield "}"

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


    result_table =  table_factory(new_data,otu_table.SampleIds,genome_table.ObservationIds,constructor=SparseGeneTable)

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
        verbose = False):
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
        verbose = False):
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
        verbose = False):
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
        curr_nsti =  float(genome_table.SampleMetadata[obs_id_idx]['NSTI'])
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
def variance_of_dot_product(variances_v1,variances_v2):
    """Calculate the variance of a dot product
    The dot product is the sum of products for each pair of items in vector1 and vector2

    The variance of each product of v1[i] and v2[i] is approximated by:
    variance_v1[i]/(v1[i]**2) + variance_v2[i]/(v2[i]**2) +2(sqrt(variance_v1[i])*sqrt(variance_v1[i]))*r(v1,v2)/(v1[i]*v2[i])
    Where r(v1,v2) is the correlation coefficient for v1 and v2.  (If v1 and v2 are uncorrelated the final term can be omitted)

    As described here: http://en.wikipedia.org/wiki/Propagation_of_uncertainty
    
    """

def variance_of_product(A,B,varA,varB,r=0):
    variance = (varA/A)**2 + (varB/B)**2 
    if r !=0: #vectors are correlated
        variance += 2*(sqrt(varA)*sqrt(varB))*r/(A*B)
    return variance

def variance_of_sum(varA,varB,r=0,sign_of_varB=1):

    variance = varA + varB + 2*(sqrt(varA)*sqrt(varB))*r*sign_of_varB
    return variance

def predict_metagenome_variances(otu_table,metagenome_table,variances):
    """Predict variances for metagenome predictions
    otu_table -- BIOM Table object of OTUs
    metagenome_table -- BIOM Table object of predicted gene counts per OTU and samples
    variances -- BIOM Table object of predicted variance in each gene count
    """

    # Each predicted gene family will be assigned a weight in the final
    # meatagenome based on the relative abundnace of the OTUs in which
    # that gene family is found.  
    metagenome_data,variance_data,overlapping_genes = extract_otu_and_genome_data(metagenome_table,\
      variances,otu_table_ids="ObservationIds",genome_table_ids="ObservationIds")
    print "OTU data:",otu_data
    print "variance_data:",variance_data
    weight_per_otu = numpy_sum(otu_data,axis=1)
    print "weight_per_otu:",weight_per_otu
    
    gene_variances = variance_of_weighted_mean(numpy_sum(array(otu_data),axis=1),array(variance_data),per_sample_axis=1)
    lower_95_CI,upper_95_CI = calc_confidence_interval_95(metagenome_table, gene_variances)
    return gene_variances, lower_95_CI,upper_95_CI
