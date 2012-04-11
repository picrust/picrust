#!/usr/bin/env python
# File created on 10 April 2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2012, The PICRUST project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 


from cogent.util.option_parsing import parse_command_line_parameters,\
  make_option
from picrust.evaluate_test_datasets import evaluate_test_dataset
from biom.parse import parse_biom_table
from picrust.parse import parse_marker_gene_copy_numbers

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","Evaluate the accuracy of all predictions in a folder, and output summary statistics.","%prog -i obs_otu_table.biom -e exp_otu_table.txt -o./evaluation_results/")]
script_info['output_description']= "Outputs will be obs,exp data points for the comparison"
script_info['required_options'] = [
 make_option('-i','--input_trait_table_fp',type="existing_filepath",help='the input trait table filepath in biom format'),
 make_option('-e','--input_expected_trait_table_fp',type="existing_filepath",help='the input expected trait table filepath'),
 make_option('-o','--output_dir',type="new_filepath",help='the output directory'),
]
script_info['optional_options'] = [
]
script_info['version'] = __version__



#def evaluate_test_dataset_dir(obs_dir_fp,exp_dir_fp):
#    """Return evaluation results from the given directories
#    
#    
#    """
#
#    #files = 
#    for f in exp_dir_fp:
#        #find obs_dir_fp
#        scatter_data_points, correlations=\
#          evaluate_test_datasets(obs_table,exp_table)



def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    if opts.verbose:
        print "Loading observed table..."
    observed_table =\
      parse_biom_table(open(opts.input_trait_table_fp,'U'))
    
    if opts.verbose:
        print "Loading expected table..."
    
    expected_table =\
      parse_biom_table(open(opts.input_expected_trait_table_fp,'U'))
    
    if opts.verbose:
        print "Calculating evaluation accuracy..."

    results = evaluate_test_dataset(observed_table,expected_table)
    print results
    #open(opts.output_otu_fp,'w').write(\
    # normalized_table.getBiomFormatJsonString('PICRUST'))


if __name__ == "__main__":
    main()
