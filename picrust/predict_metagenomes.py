#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Greg Caporaso","Jesse Zaneveld","Morgan Langille"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from numpy import abs,compress, dot, array, around, asarray,empty,zeros, sum as numpy_sum,sqrt,apply_along_axis
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
    
    # While constructing the new table, we need to  preserve metadata about the samples from the OTU table, 
    #and metadata about the gene functions from the genome table
    result_table=table_from_template(new_data,otu_table.SampleIds,genome_table.ObservationIds,\
      sample_metadata_source=otu_table,observation_metadata_source=genome_table,\
      constructor=SparseGeneTable,verbose=verbose)

    return result_table

def predict_metagenome_variances(otu_table,genome_table,\
    gene_variances,verbose=False):
    """Predict variances for metagenome predictions
    otu_table -- BIOM Table object of OTUs
    gene_table -- BIOM Table object of predicted gene counts per OTU and samples
    gene_variances -- BIOM Table object of predicted variance in each gene count
   
   Note that OTU counts are treated as constants (exactly known) rather than random variables
   for now.   If a good method for getting variance for OTU counts becomes available, this should
   be updated to treat them as random variables as well.
   """
    #Assume that OTUs are SampleIds in the genome table, but ObservationIds in the OTU table
    genome_table_otu_ids="SampleIds"
    otu_table_otu_ids="ObservationIds"
    
    #Find overlapping otus
    overlapping_otus = get_overlapping_ids(otu_table,genome_table,\
                  genome_table_ids=genome_table_otu_ids,otu_table_ids=otu_table_otu_ids)
    #Ensure they overlap fully with variance table 
    overlapping_otus = get_overlapping_ids(otu_table,gene_variances,\
                  genome_table_ids=genome_table_otu_ids,otu_table_ids=otu_table_otu_ids)
    #Filter OTU and Genome Table to contain only overlapping IDs
    #print "overlapping_otus:",overlapping_otus
    otu_table.filterObservations(lambda val,otu_id,metadata: otu_id in overlapping_otus)
    genome_table.filterSamples(lambda val,otu_id,metadata: otu_id in overlapping_otus)
    
    #Handle missing variance data
    #if gene_variances is None:
    #    gene_variances = genome_table.copy()
    #    gene_variances.transformSamples(lambda val,otu_id,metadata: val*0.0) 
    #    #TODO: test if this is faster or slower than filling numpy.zeros followed by table
    #    #construction

   
    metagenome_data = None
    metagenome_variance_data = None
    if verbose:
        print "Calculating the variance of the estimated metagenome for %i OTUs." %len(overlapping_otus)
    for otu_id in overlapping_otus:
        otu_across_samples = otu_table.observationData(otu_id)
        otu_across_genes = genome_table.sampleData(otu_id)
        otu_variance_across_genes = gene_variances.sampleData(otu_id)
        otu_contrib_to_metagenome=array([o*otu_across_genes for o in otu_across_samples])
        var_otu_contrib_to_metagenome=\
          array([scaled_variance(otu_variance_across_genes,o) for o in otu_across_samples])
        
        if metagenome_data is None:
            metagenome_data = otu_contrib_to_metagenome
            metagenome_variance_data = var_otu_contrib_to_metagenome
        else:
            metagenome_data += otu_contrib_to_metagenome
            metagenome_variance_data = variance_of_sum(metagenome_variance_data,var_otu_contrib_to_metagenome)

    data_result = metagenome_data.T    
    variance_result = metagenome_variance_data.T
    
    if verbose:
        print "Calculating metagenomic confidene intervals from variance."

    lower_95_CI,upper_95_CI=calc_confidence_interval_95(data_result,variance_result,\
      round_CI=True,min_val=0.0,max_val=None)

    
    if verbose:
        print "Generating BIOM output tables for the prediction,variance,upper confidence interval and lower confidence interval."
    
    #Wrap results into BIOM Tables
    result_data_table=\
      table_from_template(data_result,otu_table.SampleIds,\
      genome_table.ObservationIds,sample_metadata_source=otu_table,\
      observation_metadata_source=genome_table,constructor=SparseGeneTable)

    result_variance_table=\
      table_from_template(variance_result,otu_table.SampleIds,\
      genome_table.ObservationIds,sample_metadata_source=\
      otu_table,observation_metadata_source=genome_table,constructor=\
      SparseGeneTable,verbose=verbose)
    
    result_lower_CI_table=\
      table_from_template(lower_95_CI,otu_table.SampleIds,\
      genome_table.ObservationIds,sample_metadata_source=otu_table,\
      observation_metadata_source=genome_table,constructor=SparseGeneTable,\
      verbose=verbose)
    
    result_upper_CI_table=\
      table_from_template(upper_95_CI,otu_table.SampleIds,\
      genome_table.ObservationIds,sample_metadata_source=\
      otu_table,observation_metadata_source=genome_table,constructor=\
      SparseGeneTable,verbose=verbose)
    
    return result_data_table,result_variance_table,result_lower_CI_table,\
      result_upper_CI_table


def table_from_template(new_data,sample_ids,observation_ids,\
    sample_metadata_source=None,observation_metadata_source=None,\
    constructor=SparseGeneTable,verbose=False):
    """Build a new BIOM table from new_data, and transfer metadata from 1-2 existing tables"""

    #Build the BIOM table
    result_table =  table_factory(new_data,sample_ids,observation_ids,\
      constructor=SparseGeneTable)
    
    
    #Transfer sample metadata from the OTU table
    #to the metagenome table (samples are the same)
    if sample_metadata_source:
        result_table = transfer_metadata(sample_metadata_source,result_table,\
          donor_metadata_type='SampleMetadata',\
          recipient_metadata_type='SampleMetadata',verbose=verbose)
    
    #Now transfer observation metadata (e.g. gene metadata) 
    #from the genome table to the result table
    if observation_metadata_source:
        result_table = transfer_metadata(observation_metadata_source,\
          result_table,donor_metadata_type='ObservationMetadata',\
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

def variance_of_product(A,B,varA,varB,r=0):
    """Return the variance of the product of two random variables A and B"""
    variance = (varA**2)/(A**2) 
    variance += (varB**2)/(B**2) 
    #add variance due to correlation between vectors
    variance += 2*(sqrt(varA)*sqrt(varB))*r/(A*B)
    return variance

def scaled_variance(var_x,c):
    """Return the variance of a random variable X scaled by a constant c 
    
    Formula:  Var(X*c) = var(X)*(c**2)
    
    """
    return var_x*(c**2)

def variance_of_sum(varA,varB,r=0,sign_of_varB=1):

    variance = varA + varB + 2*(sqrt(varA)*sqrt(varB))*r*sign_of_varB
    return variance

   
def sum_rows_with_variance(data_array,variance_array):    
    """Sum the rows of an array, returning sum and variance arrays
        
    
    """
    if data_array.shape != variance_array.shape:
        raise ValueError("data array and variance array must have the same shape. Instead we have data.shape:%s and variance.shape:%s"%(data_array.shape,variance_array.shape))
    
    result_by_rows = None
    for row in data_array:
        if result_by_rows is None:
            result_by_rows = row
        else:
            result_by_rows += row
        last_row = row
    variance_by_rows = None
    for row in variance_array:
        if variance_by_rows is None:
            variance_by_rows = row
        else:
            new_variance = variance_of_sum(variance_by_rows,row)
            variance_by_rows = new_variance
    return result_by_rows,variance_by_rows

