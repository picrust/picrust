#!/usr/bin/env python
# File created on 09 Jan 2013
from __future__ import division

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

from cogent.util.option_parsing import parse_command_line_parameters, make_option
from biom.parse import parse_biom_table

script_info = {}
script_info['brief_description'] = "Collapse table data to a specified level in a hierarchy."
script_info['script_description'] = "This script collapses hierarchical data to a specified level. For instance, often it is useful to examine KEGG results from a higher level within the pathway hierarchy. Many genes are sometimes involved in multiple pathways, and in these circumstances (also know as a one-to-many relationship), the gene is counted for each pathway. This has a side effect of increasing the total count of genes in the table."
script_info['script_usage'] = [("","Collapse predicted metagenome results.","""%prog -i metagenome.biom -c "KEGG Pathways" -l 3 -o metagenome_at_level3.biom""")]
script_info['output_description']= "Output table is contains gene counts at a higher level within a hierarchy."
script_info['required_options'] = [\
 make_option('-i','--input_fp',type="existing_filepath",help='the predicted metagenome table'),\
 make_option('-o','--output_fp',type='new_filepath', help='the resulting table'),
 make_option('-c','--metadata_category',type='string',help='the metadata category that describes the hierarchy'),
 make_option('-l','--level',type='int',help='the level in the hierarchy to collapse to. A value of 0 is not allowed, a valuue of 1 is the highest level, and any higher value nears the leaves of the hierarchy. For instance, if the hierarchy contains 4 levels, specifying 3 would collapse at one level above being fully specified.')
]
script_info['optional_options'] = []
script_info['version'] = __version__

def make_collapse_f(category, level):
    def collapse(md):
        for path in md[category]:
            yield path[level]
    return collapse

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    if opts.level <= 0:
        parser.error("level must be greater than zero!")

    collapse_f = make_collapse_f(opts.metadata_category, opts.level)
    table = parse_biom_table(open(opts.input_fp))
    result = table.collapseObservationsByMetadata(collapse_f, one_to_many=True, 
                                                  norm=False)
    
    f = open(opts.output_fp,'w')
    f.write(result.getBiomFormatJsonString('picrust %s - categorize_by_function'\
                                           % __version__))
    f.close()

if __name__ == "__main__":
    main()
