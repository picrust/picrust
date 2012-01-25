#!/usr/bin/env python
# Author: Morgan Langille (morgan.g.i.langille@gmail.com)
# count_wagner.py

""" Application controller for Count (ancestral state resconstruction)

File created on 24 Jan 2012.

"""
from __future__ import division
import argparse
from cogent.util.table import Table
from os.path import split, splitext
from os import remove
from glob import glob
from cogent.app.util import CommandLineApplication, ResultPath, get_tmp_filename
from cogent.app.parameters import ValuedParameter, FilePath

__author__ = "Morgan Langille"
__copyright__ = "Copyright 2011-2012, The PICRUST Project"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"

class Count(CommandLineApplication):
    """ Application controller for Count an ASR tool."""

    _command = 'java -Xmx1024M -cp Count.jar ca.umontreal.iro.evolution.genecontent.AsymmetricWagner'
    _parameters = {\
     '-gain':ValuedParameter(Prefix='-',Name='gain',Delimiter=' '),\
     '-max_paralogs':ValuedParameter(Prefix='-',Name='max_paralogs',Delimiter=' ')}
    _input_handler = '_input_as_paths'
    _suppress_stdout = False
    _suppress_stderr = False

def infer_wagner_parsimony_from_paths(tree_path,trait_table_path,gain=None,max_paralogs=None,HALT_EXEC=False):
    '''Runs count application controller given path of tree and trait table and returns a Table'''
    #initialize Count app controller
    count=Count(HALT_EXEC=HALT_EXEC)

    #set the parameters
    if gain:
        count.Parameters['-gain'].on(gain)
    if max_paralogs:
        count.Parameters['-max_paralogs'].on(max_paralogs)

    #Run Count
    result = count(data=(tree_path,trait_table_path))

    #parse the results into a Cogent Table
    asr_table= parse_wagner_parsimony_output(result["StdOut"].readlines())

    return asr_table

def infer_wagner_parsimony_from_objects(tree_object,trait_table_object,gain=None,max_paralogs=None,HALT_EXEC=False):
    '''Runs count application controller given a cogent tree object and a cogent Table object'''
    pass

def parse_wagner_parsimony_output(raw_output_with_comments):
    '''Parses wagner parsimony output from Count and returns a Cogent Table object'''
    #keep only lines in with actual ASR count information
    #throw away first 2 columns and last 4 columns
    filtered_output=[x.split()[2:-4] for x in raw_output_with_comments if x[0:8] == '# FAMILY']
    #Take the first row as the header and the rest as rows in the table
    table=Table(filtered_output[0],filtered_output[1:])
    return table

def parse_command_line_parameters():
    """ Parses command line arguments """

    parser = argparse.ArgumentParser(\
        description='Uses Wagner Parsimony from Count package to infer ancestral protein counts')
    parser.add_argument('tree_file', help='a tree in newick format')
    parser.add_argument('trait_file', help='a trait table')
    parser.add_argument('--gain',dest='gain',type=float, metavar='N', help='The gain rate (default is 1)')
    parser.add_argument('--max_paralogs',dest='max_paralogs',type=int, metavar='N', help='Filter out families with more than max paralogs')
    parser.add_argument('--debug',dest='debug',default=False,action='store_true',help='Prints out command to be executed without actually running it')
    parser.add_argument('--output',dest='output',type=str,metavar='file',help='File to store output in (defaut STDOUT)')

    args = parser.parse_args()
    
    return args


if __name__ == "__main__":
    args = parse_command_line_parameters()
    
    tree_filepath = args.tree_file
    trait_table_filepath = args.trait_file
    gain = args.gain
    max_paralogs = args.max_paralogs
    
    asr_table = infer_wagner_parsimony_from_paths(\
     tree_filepath,trait_table_filepath,gain=gain,max_paralogs=max_paralogs,HALT_EXEC=args.debug)

    if args.output:
        f=open(args.output,"w")
        f.write(asr_table.tostring(sep='\t'))
        f.close()
    else:
        print asr_table.tostring(sep='\t')
    
    
