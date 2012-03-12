#!/usr/bin/env python
# File created on 22 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The PICRUST project"
__credits__ = ["Greg Caporaso","Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 
def parse_marker_gene_copy_numbers(counts_f,
                                   metadata_identifier):
    result = {}
    #skip header line
    next(counts_f)
    for line in counts_f:
        fields = line.strip().split('\t')
        try:
            #data can be floats so round them and make them integers
            copy_number = int(round(float(fields[1])))
            
        except ValueError:
            raise ValueError,\
                "Invalid type passed as copy number for OTU ID %s. Must be int-able." % (fields[0])
        if copy_number < 1:
            raise ValueError, "Copy numbers must be greater than or equal to 1."
        result[fields[0]] = {metadata_identifier:copy_number}
    return result



def parse_trait_table(trait_table_lines,delimiter="\t",has_header=True):
    """Return a header line, and a generator that will yield data fields
    
    trait_table_lines -- tab-seperated lines that may have newline characters
    
    if has_header is True, then the first non-blank, non-comment line 
    must be a header line of equal length to the number of columns, with 
    labels for the contents of each column.

    Comment lines (starting with '#') and blank lines will be ignored,
    and won't show up in the output.
    """
    header_line = None
    if not has_header:
        header_line = ''
    else:
        for i,line in enumerate(trait_table_lines):
            if not line or line.startswith("#"):
                continue
            if i == 0:
                header_line = line
                break
    if header_line is None:
        raise RuntimeError("Could not find header line in input trait table lines.  Was it skipped due to starting with a comment ('#') sign?")
    #Now that we have the header (if present) yield_trait_table_fields
    #can assume no header exists, and just parse data fields
    return header_line, yield_trait_table_fields(trait_table_lines,\
      delimiter=delimiter,has_header=False)




def yield_trait_table_fields(trait_table_lines,delimiter="\t",\
    skip_comment_lines=True,has_header=False):
    """Yield fields from trait table lines
    
    The current definition for the header lines is as follows:
    -- can't start with a comment 
    -- must be the first line
    
    """
    for i,line in enumerate(trait_table_lines):
        
        if line.startswith("#") and skip_comment_lines:
            #ignore these and remove from outputs
            continue
        
        if has_header and i == 0:
            #header line has no fields and should be skipped
            continue
        
        #Check for bad delimiters and try to intelligently warn user
        #if they used the wrong delimiter
        if delimiter not in line:
            delimiters_to_check = {"tab":"\t","space":"","comma":","}
            possible_delimiters = []
            for delim in delimiters_to_check.keys():
                if delimiters_to_check[delim] in line:
                    possible_delimiters.append(delim)
            error_line =\
                "Delimiter '%s' not in line.  The following delimiters were found:  %s.  Is the correct delimiter one of these?"
            raise RuntimeError(error_line % (delimiter,\
             ",".join(possible_delimiters)))
        

                     
        fields = line.strip().split(delimiter)
        yield fields


