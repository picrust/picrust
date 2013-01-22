#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2012, The PICRUST project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.9.0"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from numpy import dot, array, around
from biom.table import table_factory
from picrust.predict_metagenomes import get_overlapping_ids,extract_otu_and_genome_data

def partition_metagenome_contributions(otu_table,genome_table, limit_to_functions=[], remove_zero_rows=True,verbose=True):
    """Return a list of the contribution of each organism to each function, per sample
    (rewritten version using numpy)
    otu_table -- the BIOM Table object for the OTU table
    genome_table -- the BIOM Table object for the predicted genomes
    limit_to_functions -- a list of function ids to include.  If empty, include all function ids
    Output table as a list of lists with header
    Function\tOrganism\tSample\tCounts\tpercent_of_sample
    """
    
    if limit_to_functions:
        if verbose:
            print "Filtering the genome table to include only user-specified functions:",limit_to_functionss
        ok_ids = frozenset(map(str,limit_to_functions))
        
        filter_by_set = lambda vals,gene_id,metadata: str(gene_id) in ok_ids
        #filter_by_set = lambda vals,gene_id,metadata: gene_id in ok_ids
        
        #if verbose:
            #print dir(genome_table)
        #    print "Valid function ids:",genome_table.ObservationIds
        genome_table = genome_table.filterObservations(filter_by_set)
        
        if genome_table.isEmpty():
            raise ValueError("User filtering by functions (%s) removed all results from the genome table"%(str(limit_to_functions)))

    otu_data,genome_data,overlapping_ids = extract_otu_and_genome_data(otu_table,genome_table)
    #We have a list of data with abundances and gene copy numbers
    lines=[]
    result = [["Gene","Sample","OTU","GeneCountPerGenome",\
            "OTUAbundanceInSample","CountContributedByOTU",\
            "ContributionPercentOfSample","ContributionPercentOfAllSamples"]]

    #TODO refactor as array operations for speed

    #Zero-valued total counts will be set to epsilon 
    epsilon = 1e-5

    for j,gene_id in enumerate(genome_table.ObservationIds):
        all_gene_rows = []
        for k,sample_id in enumerate(otu_table.SampleIds):
            #Add raw counts for the gene in this sample to a list
            sample_gene_rows = []
            for i,otu_id in enumerate(overlapping_ids):
                otu_gene_count = genome_data[i][j]
                otu_abundance = otu_data[i][k]
                contribution =  otu_gene_count * otu_abundance
                if remove_zero_rows and contribution == 0.0:
                    #skip zero contributions
                    continue
                sample_gene_rows.append([gene_id,sample_id,otu_id,otu_gene_count,otu_abundance,contribution])
            #Now get the percentage of each genes contribution to the sample overall
            total_counts =max(epsilon,sum([float(row[-1]) for row in sample_gene_rows]))

            for row in sample_gene_rows:
                percent_of_sample = float(row[-1])/total_counts
                row.append(percent_of_sample)
            all_gene_rows.extend(sample_gene_rows)
        
        count_idx = -2 #Counts are now in the next to last position in each row
        total_counts =max(epsilon,sum([float(row[count_idx]) for row in all_gene_rows]))

        for row in all_gene_rows:
            percent_of_sample = float(row[count_idx])/total_counts
            row.append(percent_of_sample)
        lines.extend(all_gene_rows)
    result.extend(lines)

    return result

