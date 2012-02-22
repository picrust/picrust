#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The PICRUST project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 
def parse_marker_gene_copy_numbers(counts_f,
                                   metadata_identifier):
    result = {}
    for line in counts_f:
        fields = line.strip().split('\t')
        try:
            copy_number = int(fields[1])
        except ValueError:
            raise ValueError,\
             "Invalid type passed as copy number for OTU ID %s. Must be int-able." % (fields[0])
        if copy_number < 1:
            raise ValueError, "Copy numbers must be greater than or equal to 1."
        result[fields[0]] = {metadata_identifier:copy_number}
    return result