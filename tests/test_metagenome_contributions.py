#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2012, The PICRUST project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.9.1-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from biom.parse import parse_biom_table_str
from biom.table import DenseTable
from picrust.metagenome_contributions import partition_metagenome_contributions
from picrust.predict_metagenomes import predict_metagenomes,\
  calc_nsti,get_overlapping_ids, extract_otu_and_genome_data

class PartitionMetagenomeTests(TestCase):
    """ """
    
    def setUp(self):
        self.otu_table1 = parse_biom_table_str(otu_table1)
        self.genome_table1 = parse_biom_table_str(genome_table1)
        self.genome_table2 = parse_biom_table_str(genome_table2)
        self.predicted_metagenome_table1 = parse_biom_table_str(predicted_metagenome_table1)
        self.predicted_gene_partition_table = predicted_gene_partition_table
 
    def test_partition_metagenome_contributions(self):
        """partition_metagenome_contributions functions as expected with valid input"""
        #For reference, the OTU table should look like this:
        ##OTU ID Sample1 Sample2 Sample3 Sample4
        #GG_OTU_1    1.0 2.0 3.0 5.0
        #GG_OTU_2    5.0 1.0 0.0 2.0
        #GG_OTU_3    0.0 0.0 1.0 4.0

        #...and the genome table will look like this:
        ##OTU ID GG_OTU_1    GG_OTU_3    GG_OTU_2
        #f1  1.0 2.0 3.0
        #f2  0.0 1.0 0.0
        #f3  0.0 0.0 1.0
      
        #For which predict metagenomes should produce a table like this:
        ##OTU ID    Sample1 Sample2 Sample3 Sample4
        #f1  16.0    5.0 5.0 19.0
        #f2  0.0 0.0 1.0 4.0
        #f3  5.0 1.0 0.0 2.0

        #First, sanity checks
        
        #We expect to see the contributions broken down by OTU 
        metagenome_table = predict_metagenomes(self.otu_table1,self.genome_table1)
        obs = partition_metagenome_contributions(self.otu_table1,self.genome_table1)
        
        obs_text = "\n".join(["\t".join(map(str,i)) for i in obs])
        exp_text = "\n".join(["\t".join(map(str,r.split())) for r in self.predicted_gene_partition_table.split('\n')])

        #Test that the percent of all samples is always smaller than the percent of the current sample
        for l in obs[1:]:
            self.assertTrue(l[-1]<=l[-2])
        
        #Test that the summed contributions equal the metagenome table value
        sum_f1_sample1 = sum([i[5] for i in obs[1:] if (i[0]=="f1" and i[1]=="Sample1")])
        self.assertFloatEqual(sum_f1_sample1,16.0)
        sum_f2_sample1 = sum([i[5] for i in obs[1:] if (i[0]=="f2" and i[1]=="Sample1")])
        self.assertFloatEqual(sum_f2_sample1,0.0)
        sum_f3_sample1 = sum([i[5] for i in obs[1:] if (i[0]=="f3" and i[1]=="Sample1")])
        self.assertFloatEqual(sum_f3_sample1,5.0)

        for l in obs[1:]:
            gene,sample,OTU,gene_count_per_genome,otu_abundance_in_sample,count,percent,percent_all = l    
            #Test that genomes without genes don't contribute
            #Only GG_OTU_3 has f2, so for all others the gene contribution should be 0,0
            if gene == "f2" and OTU != "GG_OTU_3":
                self.assertFloatEqual(count,0.0)
                self.assertFloatEqual(percent,0.0)
                self.assertFloatEqual(percent_all,0.0)
            #Ditto for GG_OTU_2 and f3
            if gene == "f3" and OTU != "GG_OTU_2":
                self.assertFloatEqual(count,0.0)
                self.assertFloatEqual(percent,0.0)
                self.assertFloatEqual(percent_all,0.0)
            
            #Test that OTUs absent from a sample don't contribute
            if sample == "Sample1" and OTU == "GG_OTU_3":
                self.assertFloatEqual(count,0.0)
                self.assertFloatEqual(percent,0.0)
                self.assertFloatEqual(percent_all,0.0)
        
        #Having validated that this looks OK, just compare to hand-checked result
        self.assertEqual(obs_text,exp_text)
       
otu_table1 = """{"rows": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [0, 3, 5.0], [1, 0, 5.0], [1, 1, 1.0], [1, 3, 2.0], [2, 2, 1.0], [2, 3, 4.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:50:05.024661", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

genome_table1 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [1, 1, 1.0], [2, 2, 1.0]], "columns": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

genome_table2 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [1, 1, 1.0], [2, 2, 1.0]], "columns": [{"id": "GG_OTU_21", "metadata": null}, {"id": "GG_OTU_23", "metadata": null}, {"id": "GG_OTU_22", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

predicted_metagenome_table1 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 16.0], [0, 1, 5.0], [0, 2, 5.0], [0, 3, 19.0], [1, 2, 1.0], [1, 3, 4.0], [2, 0, 5.0], [2, 1, 1.0], [2, 3, 2.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T16:01:30.837052", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

predicted_gene_partition_table =\
 """Gene    Sample  OTU GeneCountPerGenome  OTUAbundanceInSample    CountContributedByOTU   ContributionPercentOfSample ContributionPercentOfAllSamples
f1  Sample1 GG_OTU_1    1.0 1.0 1.0 0.0625  0.0222222222222
f1  Sample1 GG_OTU_2    3.0 5.0 15.0    0.9375  0.333333333333
f1  Sample2 GG_OTU_1    1.0 2.0 2.0 0.4 0.0444444444444
f1  Sample2 GG_OTU_2    3.0 1.0 3.0 0.6 0.0666666666667
f1  Sample3 GG_OTU_1    1.0 3.0 3.0 0.6 0.0666666666667
f1  Sample3 GG_OTU_3    2.0 1.0 2.0 0.4 0.0444444444444
f1  Sample4 GG_OTU_1    1.0 5.0 5.0 0.263157894737  0.111111111111
f1  Sample4 GG_OTU_2    3.0 2.0 6.0 0.315789473684  0.133333333333
f1  Sample4 GG_OTU_3    2.0 4.0 8.0 0.421052631579  0.177777777778
f2  Sample3 GG_OTU_3    1.0 1.0 1.0 1.0 0.2
f2  Sample4 GG_OTU_3    1.0 4.0 4.0 1.0 0.8
f3  Sample1 GG_OTU_2    1.0 5.0 5.0 1.0 0.625
f3  Sample2 GG_OTU_2    1.0 1.0 1.0 1.0 0.125
f3  Sample4 GG_OTU_2    1.0 2.0 2.0 1.0 0.25"""

if __name__ == "__main__":
    main()
