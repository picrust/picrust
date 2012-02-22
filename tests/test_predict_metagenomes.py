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
 

from cogent.util.unit_test import TestCase, main
from biom.parse import parse_biom_table_str
from picrust.predict_metagenomes import predict_metagenomes

class PredictMetagenomeTests(TestCase):
    """ """
    
    def setUp(self):
        self.otu_table1 = parse_biom_table_str(otu_table1)
        self.genome_table1 = parse_biom_table_str(genome_table1)
        self.predicted_metagenome_table1 = parse_biom_table_str(predicted_metagenome_table1)
    
    def test_predict_metagenomes(self):
        """ predict_metagenomes functions as expected with valid input """
        actual = predict_metagenomes(self.otu_table1,self.genome_table1)
        self.assertEqual(actual,self.predicted_metagenome_table1)


otu_table1 = """{
    "id":null,
    "format": "Biological Observation Matrix v0.9",
    "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
    "type": "OTU table",
    "generated_by": "QIIME revision XYZ",
    "date": "2011-12-19T19:00:00",
    "rows":[
            {"id":"GG_OTU_1", "metadata":null},
            {"id":"GG_OTU_2", "metadata":null},
            {"id":"GG_OTU_3", "metadata":null},
            {"id":"GG_OTU_4", "metadata":null},
            {"id":"GG_OTU_5", "metadata":null}
        ],
    "columns": [
            {"id":"Sample1", "metadata":null},
            {"id":"Sample2", "metadata":null},
            {"id":"Sample3", "metadata":null},
            {"id":"Sample4", "metadata":null},
            {"id":"Sample5", "metadata":null},
            {"id":"Sample6", "metadata":null}
        ],
    "matrix_type": "sparse",
    "matrix_element_type": "int",
    "shape": [5, 6],
    "data":[[0,2,1],
            [1,0,5],
            [1,1,1],
            [1,3,2],
            [1,4,3],
            [1,5,1],
            [2,2,1],
            [2,3,4],
            [2,4,2],
            [3,0,2],
            [3,1,1],
            [3,2,1],
            [3,5,1],
            [4,1,1],
            [4,2,1]
           ]
}"""

genome_table1 = """{
    "id":null,
    "format": "Biological Observation Matrix v0.9",
    "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
    "type": "OTU table",
    "generated_by": "QIIME revision XYZ",
    "date": "2011-12-19T19:00:00",
    "rows":[
            {"id":"f1", "metadata":null},
            {"id":"f2", "metadata":null},
            {"id":"f3", "metadata":null},
            {"id":"f4", "metadata":null}
        ],
    "columns": [
            {"id":"Sample1", "metadata":null},
            {"id":"Sample2", "metadata":null},
            {"id":"Sample3", "metadata":null},
            {"id":"Sample4", "metadata":null},
            {"id":"Sample5", "metadata":null},
            {"id":"Sample6", "metadata":null}
        ],
    "matrix_type": "sparse",
    "matrix_element_type": "int",
    "shape": [4, 6],
    "data":[[0,2,1],
            [1,0,5],
            [1,1,1],
            [1,3,2],
            [1,4,3],
            [1,5,1],
            [2,2,1],
            [2,3,4],
            [2,4,2],
            [3,0,2],
            [3,1,1],
            [3,2,1],
            [3,5,1]
           ]
}"""

predicted_metagenome_table1 = """{"rows": [{"id": "f1", "metadata": null}, {"id": "f2", "metadata": null}, {"id": "f3", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 16.0], [0, 1, 5.0], [0, 2, 5.0], [0, 3, 19.0], [1, 2, 1.0], [1, 3, 4.0], [2, 0, 5.0], [2, 1, 1.0], [2, 3, 2.0]], "columns": [{"id": "Sample1", "metadata": null}, {"id": "Sample2", "metadata": null}, {"id": "Sample3", "metadata": null}, {"id": "Sample4", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2753", "matrix_type": "sparse", "shape": [3, 4], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2012-02-22T16:01:30.837052", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

if __name__ == "__main__":
    main()