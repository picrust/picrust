#!/usr/bin/env python
# Author: Morgan Langille (morgan.g.i.langille@gmail.com)
# count_wagner.py

""" Application controller for the `ace' function within the R package `ape`.

File created on Feb 2012.

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
from picrust.util import get_picrust_project_dir
from os.path import join

__author__ = "Morgan Langille"
__copyright__ = "Copyright 2011-2012, The PICRUST Project"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"


class Ace(CommandLineApplication):
    """ Application controller for 'ace' fucntion within the 'ape' R package."""

    ace_script_fp = join(get_picrust_project_dir(),'picrust','support_files','R','ace.R')
    _command = ace_script_fp
    _input_handler = '_input_as_string'
    _suppress_stdout = False
    _suppress_stderr = False


def ace_for_picrust(tree_path,trait_table_path,method='pic',HALT_EXEC=False):
    '''Runs the Ace application controller given path of tree and trait table and returns a Table'''
    #initialize Ace app controller
    ace=Ace(HALT_EXEC=HALT_EXEC)

    tmp_output_count_path=get_tmp_filename()
    tmp_output_prob_path=get_tmp_filename()
       
    as_string = " ".join([tree_path,trait_table_path,method,tmp_output_count_path,tmp_output_prob_path])
    #Run ace here
    result = ace(data=as_string)

    #Load the output into Table objects
    asr_table=LoadTable(filename=tmp_output_count_path,header=True,sep='\t')
    asr_prob_table=LoadTable(filename=tmp_output_prob_path,header=True,sep='\t')

    #Remove tmp files
    remove(tmp_output_count_path)
    remove(tmp_output_prob_path)

    return asr_table,asr_prob_table
