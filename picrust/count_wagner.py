#!/usr/bin/env python
# Author: Morgan Langille (morgan.g.i.langille@gmail.com)
# count_wagner.py

""" Application controller for Count, ancestral state resconstruction using Wagner Parsimony

File created on 24 Jan 2012.

"""
from __future__ import division
import argparse
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

class CountWagner(CommandLineApplication):
    """ Application controller for Wagner Parisomony method from the Count ASR tool."""

    _command = 'java -Xmx1024M -cp ~/Dropbox/programs/Count/Count.jar ca.umontreal.iro.evolution.genecontent.AsymmetricWagner'
    _parameters = {\
     '-gain':ValuedParameter(Prefix='-',Name='gain',Delimiter=' '),\
     '-max_paralogs':ValuedParameter(Prefix='-',Name='max_paralogs',Delimiter=' ')}
    _input_handler = '_input_as_paths'
    _suppress_stdout = False
    _suppress_stderr = False

def infer_ancestral_counts_from_paths(tree_path,trait_table_path,gain=None,max_paralogs=None,HALT_EXEC=False):

    count=CountWagner(HALT_EXEC=HALT_EXEC)
    if gain:
        count.Parameters['-gain'].on(gain)

    if max_paralogs:
        count.Parameters['-max_paralogs'].on(max_paralogs)
    
    result = count(data=(tree_path,trait_table_path))
    
    print "StdOut:",result["StdOut"].read()
    print "StdErr:",result["StdErr"].read()
    print "Return code:",result["ExitStatus"] 
    
    #results = parse_reconstruction_output(bayestraits_result['StdOut'].readlines())
    #print "Reconstructions:",results
    #outfile='count_results.txt'
    #Reconstruction results
    #f = open(outfile,"w+")
    #f.writelines(results)
    #f.close()
    return ''

def infer_ancestral_counts_from_objects(tree_object,trait_table_object,gain=None,max_paralogs=None,HALT_EXEC=False):
    pass

def parse_command_line_parameters():
    """ Parses command line arguments """

    parser = argparse.ArgumentParser(\
        description='Uses Wagner Parsimony from Count package to infer ancestral protein counts')
    parser.add_argument('tree_file', help='a tree in newick format')
    parser.add_argument('trait_file', help='a trait table')
    parser.add_argument('--gain',dest='gain',type=float, metavar='N', help='The gain rate (default is 1)')
    parser.add_argument('--max_paralogs',dest='max_paralogs',type=int, metavar='N', help='Filter out families with more than max paralogs')
    parser.add_argument('--debug',dest='debug',default=False,action='store_true',help='Prints out command to be executed without actually running it')

    args = parser.parse_args()
    
    return args


if __name__ == "__main__":
    args = parse_command_line_parameters()
    
    tree_filepath = args.tree_file
    trait_table_filepath = args.trait_file
    gain = args.gain
    max_paralogs = args.max_paralogs
    
    asr_trait_table_path = infer_ancestral_counts_from_paths(\
     tree_filepath,trait_table_filepath,gain=gain,max_paralogs=max_paralogs,HALT_EXEC=args.debug)
    
    
