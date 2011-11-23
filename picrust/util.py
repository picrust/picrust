#!/usr/bin/env python
# File created on 23 Nov 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The PICRUST project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from os.path import abspath, dirname

def get_picrust_project_dir():
    """ Returns the top-level PICRUST directory
    """
    # Get the full path of util.py
    current_file_path = abspath(__file__)
    # Get the directory containing util.py
    current_dir_path = dirname(current_file_path)
    # Return the directory containing the directory containing util.py
    return dirname(current_dir_path)
