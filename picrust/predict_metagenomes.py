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



def dot_product_with_variance(v1,v2,variances_v1,variances_v2,epsilon=1e-20):
    """Calculate the variance of a dot product
    The dot product is the sum of products for each pair of items in vector1 and vector2
    """
    #For now assume 2D matrices!
    if v1.ndim != 2 or v2.ndim !=2:
        raise ValueError("Currently only multiplication of 2D matrices is supported")
    #Ensure data matrices are the same shape
    if v1.shape[0] != v2.shape[1]:
        raise ValueError("Matrix multiplication requires an equal number of rows in matrix A as columns in matrix B. A:%s. B:%s"%(v1.shape,v2.shape))
    #Ensure variance matrices have the same shape
    #if  variances_v1.shape[0] != variances_v2.shape[1]:
    #    raise ValueError("Matrix multiplication with variance requires an equal number of rows in variance matrix A as columns in variance matrix B. A:%s. B:%s"%(variances_v1.shape,variances_v2.shape))
    #Ensure variance and data matrix have the same shape
    #if variances_v1.shape[0] != v2.shape[1]:
    #    raise ValueError("Matrix multiplication with variance requires an equal number of rows in variance matrix A as columns in data matrix B. A:%s. B:%s"%(variances_v1.shape,v2.shape))
    
    print "v1:",v1
    print "v1.T:",v1.T
    print "v2:",v2
    #v1 = v1.T
    data_result= v1.dot(v2)

    #First get the dot product result.
    #Each entry in the final dot product simply results from multiplication and summing of 
    #elements in v1 and v2
    
    #X-axis contribution 
    print "final shape:",data_result.shape
    print "v1.ravel",v1.ravel()
    print "v2.ravel",v2.ravel()
    
    result = sum([v1[i,:]*v2[:,i] for i in range(v1.shape[0])])
    print "simple result:",result
    variance_result = None
    return data_result,variance_result




def variance_of_vector_product(v1,v2,var_v1,var_v2,r=0):
    """calculate the variance of the product of two 1D vectors"""
    if v1.ndim != 1 or v2.ndim !=1:
        raise ValueError("variance of vector product can only accept two 1D vectors")
    if v1.shape != v2.shape or var_v1.shape != var_v2.shape or var_v1.shape != v1.shape:
        raise ValueError("Vectors must be the same shape")
    #multiply elements
    return variance_of_product(v1,v2,var_v1,var_v2,r)

def variance_of_product(A,B,varA,varB,r=0):
    variance = (varA/A)**2 + (varB/B)**2 
    #add variance due to correlation between vectors
    variance += 2*(sqrt(varA)*sqrt(varB))*r/(A*B)
    return variance

def variance_of_sum(varA,varB,r=0,sign_of_varB=1):

    variance = varA + varB + 2*(sqrt(varA)*sqrt(varB))*r*sign_of_varB
    return variance


def predict_metagenome_variances(otu_table,genome_table,\
    gene_variances=None,otu_variances=None,epsilon=1e-20):
    """Predict variances for metagenome predictions
    otu_table -- BIOM Table object of OTUs
    gene_table -- BIOM Table object of predicted gene counts per OTU and samples
    gene_variances -- None or BIOM Table object of predicted variance in each gene count
      If set to None, an empty array containing epsilon (very small) variance will be created 
      in the same shape as the gene table data.
    otu_variances -- None or BIOM Table object of predicted variance in each OTU count
      If set to None, an empty array containing epsilon (very small) variance will be created 
      in the same shape as the otu table data.
   """
    #Assume that OTUs are SampleIds in the genome table, but ObservationIds in the OTU table
    genome_table_otu_ids="SampleIds"
    otu_table_otu_ids="ObservationIds"
    
    #Find overlapping otus
    overlapping_otus = get_overlapping_ids(otu_table,genome_table,\
                  genome_table_ids=genome_table_otu_ids,otu_table_ids=otu_table_otu_ids)
   
    #Filter OTU and Genome Table to contain only overlapping IDs
    print "overlapping_otus:",overlapping_otus
    otu_table.filterObservations(lambda val,otu_id,metadata: otu_id in overlapping_otus)
    genome_table.filterSamples(lambda val,otu_id,metadata: otu_id in overlapping_otus)
    
    #Handle missing variance data
    if gene_variances is None:
        gene_variances = genome_table.copy()
        gene_variances.transformSamples(lambda val,otu_id,metadata: val*0.0) 
        #TODO: test if this is faster or slower than filling numpy.zeros followed by table
        #construction

    if otu_variances is None:
        otu_variances = otu_table.copy()
        otu_variances.transformObservations(lambda val,otu_id,metadata: val*0.0) #just setting to zero doesn't work! 
        #TODO: test if this is faster or slower than filling numpy.zeros followed by table
        #construction


    metagenome_data = None
    metagenome_variance_data = None
    for otu_id in overlapping_otus:
        otu_across_samples = otu_table.observationData(otu_id)
        print "otu_across_samples:",otu_across_samples
        otu_across_genes = genome_table.sampleData(otu_id)
        print "otu_across_genes:",otu_across_genes
        
        otu_variance_across_samples = otu_variances.observationData(otu_id)
        print "var(otu_across_samples):",otu_variance_across_samples
        otu_variance_across_genes = gene_variances.sampleData(otu_id)
        print "var(otu_across_genes):",otu_variance_across_genes
        otu_contrib_to_metagenome=array([o*otu_across_genes for o in otu_across_samples])
        print "OTU contribution to metagenome:",otu_contrib_to_metagenome
        S = otu_across_samples+epsilon
        G = otu_across_genes+epsilon
        var_otu_contrib_to_metagenome=\
          array([variance_of_product(o,G,otu_variance_across_samples[i],otu_variance_across_genes,r=0.0)\
          for i,o in enumerate(S)])
        print "var(OTU contribution to metagenome:)",var_otu_contrib_to_metagenome
        if metagenome_data is None:
            metagenome_data = otu_contrib_to_metagenome
            metagenome_variance_data = var_otu_contrib_to_metagenome
        else:
            metagenome_data += otu_contrib_to_metagenome
            metagenome_variance_data = variance_of_sum(metagenome_variance_data,var_otu_contrib_to_metagenome)

    print metagenome_data.T
    print predict_metagenomes(otu_table,genome_table)
    data_result = metagenome_data.T    
    variance_result = metagenome_variance_data.T
    print variance_result
    #TODO: wrap into BIOM Table
    return data_result,variance_result

def predict_metagenome_variances_backup(otu_table,genome_table,\
    gene_variances=None,otu_variances=None,epsilon=1e-20):
    """Predict variances for metagenome predictions
    otu_table -- BIOM Table object of OTUs
    gene_table -- BIOM Table object of predicted gene counts per OTU and samples
    gene_variances -- None or BIOM Table object of predicted variance in each gene count
      If set to None, an empty array containing epsilon (very small) variance will be created 
      in the same shape as the gene table data.
    otu_variances -- None or BIOM Table object of predicted variance in each OTU count
      If set to None, an empty array containing epsilon (very small) variance will be created 
      in the same shape as the otu table data.

    Method description:

    Metagenome prediction is calculated using the dot product
    of the transposed OTU table and the genome table
    
    To calculate the effects of this operation on variance in the gene product,
    we need to find the intermediate steps, and account for the variance steps of 
    each.

    From the  numpy 1.7 documentation(dot function, http://docs.scipy.org/doc/numpy/reference/generated/numpy.dot.html):
    'For 2-D arrays it is equivalent to matrix multiplication, and for 1-D arrays to inner product of vectors (without complex conjugation). For N dimensions it is a sum product over the last axis of a and the second-to-last of b:

    dot(a, b)[i,j,k,m] = sum(a[i,j,:] * b[k,:,m])'

    Therefore, the strategy employed is to decompose the dot product step in predict_metagenomes.py into smaller steps
    that each consist of only addition or matrix multiplication for 1D arrays.  Uncertainty is then propagated across
    each of these smaller steps.
    """
     
    #Restrict to informative  ids 
    otu_data,genome_data,overlapping_otus =\
      extract_otu_and_genome_data(otu_table,genome_table)
    
    otu_data = array(otu_data).T
    #otu_data=array(otu_data)
    genome_data = array(genome_data)
    #genome_data = array(genome_data).T
    #Create a default array with tiny variance if none is supplied
    
    if otu_variances is None:
        otu_variance_data = zeros(otu_data.shape)
        #otu_variance_data.fill(epsilon)
    else:
        otu_variance_data = array(otu_variances)
    
    if gene_variances is None:
        gene_variance_data = zeros(genome_data.shape)
        #gene_Variance_data.fill(epsilon)
    else:
        #We only care about the new variance data
        redundant_otu_data,gene_variance_data,redundant_overlapping_otus =\
          extract_otu_and_genome_data(otu_table,gene_variances)
        gene_variance_data = array(gene_variance_data)
   
    print "otu_data.shape:",otu_data.shape
    print "genome_data.shape:",genome_data.shape
    results,variances = dot_product_with_variance(otu_data,\
      genome_data,otu_variance_data,gene_variance_data,epsilon=epsilon)
    
    #print "varaince_result:",variance_result
    #If we are rounding results, round them:
    data_result=around(results.T)
    variance_result = around(variances.T)
    
    #Convert to BIOM Table output
    #print "Genome table, sampleIDs:",genome_table.SampleIds
    #print "Genome table, obsIDs:",genome_table.ObservationIds
    #print "gene result:",data_result.shape
    #print "variance result:",variance_result.shape

    data_result_table =  table_factory(data_result,otu_table.SampleIds,genome_table.SampleIds,constructor=SparseGeneTable)
    variance_result_table =  table_factory(variance_result,otu_table.SampleIds,genome_table.SampleIds,constructor=SparseGeneTable)
    return data_result_table, variance_result_table
    
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
    print "result_by_rows:",result_by_rows 
    variance_by_rows = None
    for row in variance_array:
        print "row in variance array:",row 
        if variance_by_rows is None:
            variance_by_rows = row
        else:
            print "pre-sum variance_by_rows:",variance_by_rows
            new_variance = variance_of_sum(variance_by_rows,row)
            print "new_variance:",new_variance
            variance_by_rows = new_variance
            print "variance_by_rows:",variance_by_rows
    return result_by_rows,variance_by_rows

