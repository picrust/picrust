#!/usr/bin/env python
# Author: Morgan Langille (morgan.g.i.langille@gmail.com)
# count_wagner.py

""" Application controller for the `ace' function within the R package `ape`.

File created on Feb 2012.

"""
from __future__ import division
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
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Morgan Langille", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.9.0-dev"
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

    # Overridden to call script with R rather than directly - this is useful
    # because permisssions on the script are set to 644 when PICRUSt is installed
    # with setup.py. This is fine if we're executing it with R, but not if we're
    # trying to execute it directly.
    def _get_base_command(self):
        """ Returns the full command string 

            input_arg: the argument to the command which represents the input 
                to the program, this will be a string, either 
                representing input or a filename to get input from
         """
        command_parts = []
        # Append a change directory to the beginning of the command to change 
        # to self.WorkingDir before running the command
        # WorkingDir should be in quotes -- filenames might contain spaces
        cd_command = ''.join(['cd ',str(self.WorkingDir),';'])
        if self._command is None:
            raise ApplicationError, '_command has not been set.'
        command = self._command
        parameters = self.Parameters
        
        command_parts.append(cd_command)
        command_parts.append("R")
        command_parts.append("-f")
        command_parts.append(command)
        command_parts.append("--args")
        command_parts.append(self._command_delimiter.join(filter(\
            None,(map(str,parameters.values())))))
      
        return self._command_delimiter.join(command_parts).strip()
    BaseCommand = property(_get_base_command)

def ace_for_picrust(tree_path,trait_table_path,method='pic',HALT_EXEC=False):
    '''Runs the Ace application controller given path of tree and trait table and returns a Table'''
    #initialize Ace app controller
    ace=Ace(HALT_EXEC=HALT_EXEC)

    tmp_output_count_path=get_tmp_filename()
    tmp_output_prob_path=get_tmp_filename()

    #quote file names
    tree_path='"{0}"'.format(tree_path)
    trait_table_path='"{0}"'.format(trait_table_path)
    
    as_string = " ".join([tree_path,trait_table_path,method,tmp_output_count_path,tmp_output_prob_path])
    #Run ace here
    result = ace(data=as_string)

    #Load the output into Table objects
    try:
        asr_table=LoadTable(filename=tmp_output_count_path,header=True,sep='\t')
    except IOError:
        raise RuntimeError,\
         ("R reported an error on stderr:"
          " %s" % "\n".join(result["StdErr"].readlines()))
    
    asr_prob_table=LoadTable(filename=tmp_output_prob_path,header=True,sep='\t')

    #Remove tmp files
    remove(tmp_output_count_path)
    remove(tmp_output_prob_path)

    return asr_table,asr_prob_table
