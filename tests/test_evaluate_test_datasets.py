#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2012, The PICRUST project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from picrust.evaluate_test_datasets import convert_vals_to_spearman_ranks,\
spearman_correlation, calc_spearman_t,evaluate_test_dataset

from biom.parse import parse_biom_table_str

class EvaluateTestDatasetTests(TestCase):
    """ """
    
    def setUp(self):
       self.genome_table1 = parse_biom_table_str(genome_table1)
       self.genome_table2 = parse_biom_table_str(genome_table2)
    
    def test_evaluate_test_datasets(self):
        """evalutate_test_datasets returns data points and correlations"""
        print self.genome_table1.ObservationIds
        print self.genome_table2.ObservationIds
        
        obs= evaluate_test_dataset(self.genome_table1,self.genome_table2)
        print "Obs:",obs


    def test_convert_vals_to_spearman_ranks(self):
        """convert_vals_to_spearman_ranks converts a list of floats to their relative ranks."""
        
        #Example from Spearman Wikipedia page
        #http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient
        #TODO: add some examples from more formal sources
        
        
        unordered_vals = [1.2,0.8,1.2,18,2.3]
        exp_ranks = [3.5,5,3.5,1,2]
        
        obs = convert_vals_to_spearman_ranks(unordered_vals)
        self.assertFloatEqual(obs,exp_ranks)

    def test_spearman_correlation(self):
        """spearman_correlation calculates Spearman rank correlations"""
        
        #Test data taken from Wikipedia:
        #http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient
        
        x_vals = [106,86,100,101,99,103,97,113,112,110]
        y_vals = [7,0,27,50,28,29,20,12,6,17]
        
        r,prob = spearman_correlation(x_vals,y_vals,tails='high')  
        exp_rho = -0.175757575
        exp_prob_t_high = 0.6864058
        self.assertFloatEqual(r,exp_rho)
        self.assertFloatEqual(prob,exp_prob_t_high)

    def test_calc_spearman_t(self):
        """calc_spearman_t should produce an adjusted t statistic"""
        r = -0.175757575
        n = 10
        exp_prob_t_high = 0.6864058
        obs_prob_t_high = calc_spearman_t(r,n,'high')
        
        self.assertFloatEqual(obs_prob_t_high,exp_prob_t_high)

genome_table1 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [1, 1, 1.0], [2, 2, 1.0]], "columns": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

genome_table2 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 3.0], [0, 2, 4.0], [1, 1, 1.0], [2, 2, 2.0]], "columns": [{"id": "GG_OTU_21", "metadata": null}, {"id": "GG_OTU_23", "metadata": null}, {"id": "GG_OTU_22", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

if __name__ == "__main__":
    main()
