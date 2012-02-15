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
from os import remove, environ
from glob import glob
from cogent.app.util import CommandLineApplication, ResultPath, get_tmp_filename
from cogent.app.parameters import ValuedParameter, FilePath
from cogent import LoadTree
from cogent import LoadTable

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
     
    count_fp = environ['COUNT_JAR']
    _command = 'java -Xmx1024M -cp ' + count_fp + ' ca.umontreal.iro.evolution.genecontent.AsymmetricWagner'
    _parameters = {\
     '-gain':ValuedParameter(Prefix='-',Name='gain',Delimiter=' '),\
     '-max_paralogs':ValuedParameter(Prefix='-',Name='max_paralogs',Delimiter=' ')}
    _input_handler = '_input_as_paths'
    _suppress_stdout = False
    _suppress_stderr = False

    #need to overide this method since using command that is not executable
    def _error_on_missing_application(self,params):
        pass

def wagner_for_picrust(tree_path,trait_table_path,gain=None,max_paralogs=None,HALT_EXEC=False):
    '''Runs count application controller given path of tree and trait table and returns a Table'''
    #initialize Count app controller
    count=Count(HALT_EXEC=HALT_EXEC)

    #set the parameters
    if gain:
        count.Parameters['-gain'].on(gain)
    if max_paralogs:
        count.Parameters['-max_paralogs'].on(max_paralogs)

    #Create a transposed trait table and store it as a tmp file
    table = LoadTable(filename=trait_table_path,header=True,sep='\t')
    table = table.transposed(new_column_name=table.Header[0])
    tmp_table_path =get_tmp_filename()
    table.writeToFile(tmp_table_path,sep='\t')
       
    #Run Count here
    result = count(data=(tree_path,tmp_table_path))

    #Remove tmp file
    remove(tmp_table_path)

    tree=LoadTree(tree_path)
    #get the number of tips in the tree
    num_of_tips=sum(1 for x in tree.iterTips())

    #parse the results into a Cogent Table
    asr_table= parse_wagner_parsimony_output(result["StdOut"].readlines(),remove_num_tips=num_of_tips)

    #transpose the table
    asr_table = asr_table.transposed(new_column_name=asr_table.Header[0])

    return asr_table

def infer_wagner_parsimony_from_objects(tree_object,trait_table_object,gain=None,max_paralogs=None,HALT_EXEC=False):
    '''Runs count application controller given a cogent tree object and a cogent Table object'''
    pass

def parse_wagner_parsimony_output(raw_output_with_comments,remove_num_tips=0):
    '''Parses wagner parsimony output from Count and returns a Cogent Table object'''

    #keep only lines with actual ASR count information
    #throw away first 2 columns and last 4 columns (these are extra stuff from Count)
    filtered_output=[x.split()[2:-4] for x in raw_output_with_comments if x[0:8] == '# FAMILY']

    if(remove_num_tips):
        #remove columns that contain trait data for tips (not internal node data) 
        filtered_output=[[x[0]]+ x[remove_num_tips+1:] for x in filtered_output]
    
    #Take the first row as the header and the rest as rows in the table
    table=Table(filtered_output[0],filtered_output[1:])

    return table
