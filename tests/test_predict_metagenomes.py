#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2012, The PICRUST project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from biom.parse import parse_biom_table_str
from biom.table import DenseTable
from picrust.predict_metagenomes import predict_metagenomes,\
  calc_nsti,get_overlapping_ids,\
  extract_otu_and_genome_data

class PredictMetagenomeTests(TestCase):
    """ """
    
    def setUp(self):
        self.otu_table1 = parse_biom_table_str(otu_table1)
        self.genome_table1 = parse_biom_table_str(genome_table1)
        self.genome_table2 = parse_biom_table_str(genome_table2)
        self.predicted_metagenome_table1 = parse_biom_table_str(predicted_metagenome_table1)
 
    def test_predict_metagenomes(self):
        """ predict_metagenomes functions as expected with valid input """
        actual = predict_metagenomes(self.otu_table1,self.genome_table1)
        self.assertEqual(actual,self.predicted_metagenome_table1)

    def test_predict_metagenomes_value_error(self):
        """ predict_metagenomes raises ValueError when no overlapping otu ids """
        self.assertRaises(ValueError,predict_metagenomes,self.otu_table1,self.genome_table2)


      
otu_table1 = """{"rows": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [0, 3, 5.0], [1, 0, 5.0], [1, 1, 1.0], [1, 3, 2.0], [2, 2, 1.0], [2, 3, 4.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:50:05.024661", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

genome_table1 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [1, 1, 1.0], [2, 2, 1.0]], "columns": [{"id": "GG_OTU_1", "metadata": null}, {"id": "GG_OTU_3", "metadata": null}, {"id": "GG_OTU_2", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

genome_table2 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [1, 1, 1.0], [2, 2, 1.0]], "columns": [{"id": "GG_OTU_21", "metadata": null}, {"id": "GG_OTU_23", "metadata": null}, {"id": "GG_OTU_22", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T20:49:58.258296", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

predicted_metagenome_table1 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 16.0], [0, 1, 5.0], [0, 2, 5.0], [0, 3, 19.0], [1, 2, 1.0], [1, 3, 4.0], [2, 0, 5.0], [2, 1, 1.0], [2, 3, 2.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T16:01:30.837052", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

if __name__ == "__main__":
    main()
