#!/usr/bin/env python

"""name_mapping.py


This library contains functions for mapping strain names provided
by 3rd party resources onto NCBI taxon IDs


Matches are performed by:
    a. The NCBI name
    b. The aliases of the NCBI name

"""

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2007-2011, The PICRUST Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

import os
from sys import argv
from cogent.parse.ncbi_taxonomy import NcbiTaxonomyFromFiles,\
  NcbiNameParser,NcbiTaxonParser,NcbiTaxonLookup

def loose_ncbi_name_lookup(names,strict=False):
    """Returns dict mapping taxon id -> NCBI scientific name."""
    result = {}
    for name in names:
        if name.NameClass in ['scientific name',\
          'equivalent name','synonym','genbank synonym','misspelling',\
          'genbank common name']:
            result[name.Name] = (name.UniqueName or name.Name,name.TaxonId)
            result[name.UniqueName] = (name.UniqueName or name.Name,name.TaxonId)
    return result

def loose_ncbi_name_lookup_from_file(names_file):
    """Returns name dictionary from names files."""
    names = loose_ncbi_name_lookup(NcbiNameParser(names_file))
    return names
                                                

def exact_match_name(name_dict,genus,species,strain):
    """Look for an exact match of a name in the NCBI taxonomy"""
    name = " ".join([genus,species,strain])
    ncbi_name = name_dict.get(name,None)
    return ncbi_name

def loose_match_name(name_dict,genus,species,strain,\
        handle_rrndb_names = True):
    name_pieces = [genus,species,strain]
    
    
    # handle rrndb specific weird strain designation
    if '"' in strain and handle_rrndb_names:
        fixed, variable,fixed2 =strain.split('"',2)
        variable_options = variable.split(",")
        for option in variable_options:
            strain = "".join([fixed,option,fixed2])
            print "trying strain option:",strain
            name =  loose_match_name(name_dict,\
                genus,species,strain)
            print "Got:", name
            if name is not None:
                return name
    
    if '(' in strain and handle_rrndb_names:
        start,middle,stop = strain.split(["(",")"])
        

    exact_match = exact_match_name(name_dict,genus,species,strain)
    if exact_match is not None:
        return exact_match
    
    for i in range(len(name_pieces)-1):
            partial_name = " ".join(name_pieces[:-i])
            if not partial_name.strip():
                return None
            
            ncbi_name = name_dict.get(partial_name.get,None)
            if ncbi_name is not None:
                return ncbi_name

            #If not try substitutions


def main(argv):
    """Generate test trees given parameters"""
    file_with_names = open(argv[1],"U")  # temp 
    
    non_matched_names = []
    matched_names = []
    name_dict = loose_ncbi_name_lookup_from_file(open("./names.dmp"))
    #rrndb example
    for line in file_with_names:
        fields = line.split("\t")
        genus = fields[0]
        species = fields[1]
        strain = fields[2]
        # handle rrndb specific weird strain designation
    
        ncbi_name = loose_match_name(name_dict,genus,species,strain) 
        if ncbi_name is not None:
            matched_names.append(ncbi_name)
        else:
            non_matched_names.append(" ".join([genus,species,strain]))

    print len(matched_names)
    print len(non_matched_names)
    print non_matched_names

if __name__ == "__main__":
    main(argv)
