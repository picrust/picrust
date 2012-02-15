#!/usr/bin/env python
# File created on 15 Jul 2011
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011, The PICRUST Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from os.path import splitext
from cogent import LoadTree
from cogent.util.option_parsing import parse_command_line_parameters,\
    make_option



def reformat_tree_and_trait_table(tree,trait_table_lines,\
    input_trait_table_delimiter="\t", output_trait_table_delimiter="\t",\
    filter_table_by_tree_tips=True, convert_trait_floats_to_ints=False,\
    filter_tree_by_table_entries=True,convert_to_bifurcating=False,\
    add_branch_length_to_root=False, name_unnamed_nodes=True,min_branch_length=0.0001):
    """Return a full reformatted tree,pruned reformatted tree  and set of trait table lines 

    tree - a PyCogent PhyloNode tree object
    trait_table_lines -- the lines of a trait table, where 
      the rows are organisms and the columns are traits (e.g. gene counts).

    This function combines the various reformatting functions in the 
    library into a catch-all reformatter.  
    """
    
    

    input_tree = tree
    #Name unnamed nodes
    if name_unnamed_nodes:
        input_tree.nameUnnamedNodes()
    
    
    #map trait table ids to tree ids
    

    #Then filter to include only
    if filter_table_by_tree_tips:
        trait_table_lines = filter_table_by_presence_in_tree(input_tree,\
          trait_table_lines,delimiter=input_trait_table_delimiter)

    #Convert floating point values to ints
    if convert_trait_floats_to_ints:
        trait_table_lines = convert_trait_values(\
          trait_table_lines,conversion_fn = int,delimiter=input_trait_table_delimiter)
   
    #Write out results
    #output_trait_table_file.writelines(trait_table_lines)
    #trait_table.close()
    #output_trait_table_file.close()
   
    # Use the reformatted trait table to filter tree
    #trait_table = open(output_table_fp,"U")
    #trait_table_lines = trait_table.readlines()

    if filter_tree_by_table_entries:
        input_tree = filter_tree_tips_by_presence_in_table(input_tree,\
          trait_table_lines,delimiter=input_trait_table_delimiter)

    # Tree reformatting

    if convert_to_bifurcating:
        input_tree = input_tree.bifurcating() # Required by most ancSR programs


    if convert_to_bifurcating:
        input_tree = ensure_root_is_bifurcating(input_tree)
        # The below nutty-looking re-filtering step is necessary
        # When ensuring the root is bifurcating, internal nodes can get moved to the tips
        # So without additional filtering we get unannotated tip nodes
        if filter_tree_by_table_entries:
            input_tree = filter_tree_tips_by_presence_in_table(input_tree,\
              trait_table_lines,delimiter=input_trait_table_delimiter)

    if min_branch_length:
        input_tree = set_min_branch_length(input_tree,min_length = min_branch_length)

    if add_branch_length_to_root:
        input_tree = add_branch_length_to_root(input_tree,root_name=input_tree.Name,\
          root_length=min_branch_length)

    input_tree.prune()

    return input_tree, trait_table_lines


def nexus_lines_from_tree(tree):
    """Return NEXUS formatted lines from a PyCogent PhyloNode tree"""
    lines = ["#NEXUS"]
    lines.extend(make_nexus_trees_block(tree))
    return lines

def add_branch_length_to_root(tree, root_name ="root",root_length=0.0001):
    """Add branch length to the root of a tree if it's shorter than root_length
    tree -- A PyCogent PhyloNode object
    root_name -- the name of the root node
    root_length -- the desired minimum root length
    This is required by some programs such as BayesTraits"""
    
    root = tree.getNodeMatchingName(root_name)
    root.Length = max(root.Length,root_length)
    return tree 
        

def set_min_branch_length(tree,min_length= 0.0001):
    """Return tree modified so that all branchlengths are >= min_length.

    tree -- a PyCogent PhyloNode object"""

    for node in tree.preorder():
        if not node.Parent:
            continue
        node.Length = max(node.Length,min_length)
    return tree


def make_nexus_trees_block(tree):
    """Generate a NEXUS format 'trees' block for a given tree
    
    WARNING:  Removes names from internal nodes, as these cause problems
    downstream
    """

    # First generate the mappings for the NEXUS translate command
    trees_block_template =\
      ["begin trees;",\
      "\ttranslate"]
    name_mappings = {}
    line = None
    for i,node in enumerate(tree.iterTips()):
        name_mappings[node.Name] = i
        if line:
            trees_block_template.append(line)
        
        line = "\t\t%i %s," %(i,node.Name)
    # The last line needs a semicolon rather than a comma
    line = "\t\t%i %s;" %(i,node.Name)
    trees_block_template.append(line)
    
    
    # Reformat tree newick such that names match NEXUS translation table
    for name_to_fix in name_mappings.keys():
        node_to_rename = tree.getNodeMatchingName(name_to_fix)
        node_to_rename.Name=name_mappings[name_to_fix]
    for nonTipNode in tree.iterNontips():
        nonTipNode.Name=''
    

    
    tree_newick = tree.getNewick(with_distances=True)
    #for name_to_fix in name_mappings.keys():
    #    tree_newick = tree_newick.replace(name_to_fix+":",str(name_mappings[name_to_fix])+":")
    #for nonTipNode in tree.iterNontips():
    #    tree_newick = tree_newick.replace(nonTipNode.Name+":","")
    #tree_newick = tree_newick.replace(root_name,"")


    tree_template  = "\t\ttree %s = %s" # tree name then newick string
    line = tree_template % ("PyCogent_tree",tree_newick)
    trees_block_template.append(line)
    
    trees_block_template.append("end;")
    return trees_block_template


def filter_table_by_presence_in_tree(tree,trait_table_lines,name_field_index = 0,delimiter="\t"):
    """yield lines of a trait table lacking organisms missing from the tree"""
    tree_tips = [node.Name.strip() for node in tree.preorder()]
    
    for fields in yield_trait_table_fields(trait_table_lines,delimiter):
        curr_name = fields[name_field_index].strip()
        if curr_name not in tree_tips:
            continue
        yield delimiter.join(fields)+"\n"

def convert_trait_values(trait_table_lines,name_field_index=0,delimiter="\t",conversion_fn = int):
    """Convert trait values by running conversion_fn on each"""
    for fields in yield_trait_table_fields(trait_table_lines,delimiter): 
        new_fields = []
        for i,field in enumerate(fields):
            if i != name_field_index:
                new_fields.append(str(conversion_fn(float(field))))
            else:
                new_fields.append(field)
        yield delimiter.join(new_fields)+"\n"





def yield_trait_table_fields(trait_table_lines,delimiter="\t",skip_comment_lines=True,max_field_len=100):
    """Yield fields from trait table lines"""
    for line in trait_table_lines:
        #print "Parsing line:\n",line[0:min(100,len(line))],"..."
        if line.startswith("#") and skip_comment_lines:
            continue
        
        if delimiter not in line:
            delimiters_to_check = {"tab":"\t","space":"","comma":","}
            possible_delimiters = []
            for delim in delimiters_to_check.keys():
                if delimiters_to_check[delim] in line:
                    possible_delimiters.append(delim)
            error_line = "Delimiter '%s' not in line.  The following delimiters were found:  %s.  Is the correct delimiter one of these?"
            raise RuntimeError(error_line % (delimiter,",".join(possible_delimiters)))
                     
        fields = line.split(delimiter)
        yield fields



def ensure_root_is_bifurcating(tree,root_name='root'):
    """Remove child node of root if it is a single child"""
    root_node = tree.getNodeMatchingName(root_name)
    if len(root_node.Children) == 1:
        print "rerooting to avoid monotomy at root"
        tree = tree.rootedAt(root_node.Children[0].Name)
        #tree.remove(root_node)
    tree.prune()

    return tree

def filter_tree_tips_by_presence_in_table(tree,trait_table_lines,name_field_index = 0,delimiter="\t"):
    """yield a tree lacking organisms missing from the trait table"""
    org_ids_in_trait_table = []
    new_tree = tree.deepcopy()
    
    for fields in yield_trait_table_fields(trait_table_lines, delimiter):
        curr_org = fields[name_field_index].strip()
        org_ids_in_trait_table.append(curr_org)
    

    # Build up a list of tips to prune
    tips_to_prune = []
    n_tips_not_to_prune = 0
    for tip in tree.iterTips():
        if tip.Name.strip() not in org_ids_in_trait_table:
            tips_to_prune.append(tip.Name)
            #print tip.Name
            #print org_ids_in_trait_table[0]
        else:
            n_tips_not_to_prune += 1
    if not n_tips_not_to_prune:
        raise RuntimeError(\
          "filter_tree_tips_by_presence_in_table:  operation would remove all tips.  Is this due to a formatting error in inputs?")
    
    # print "Tips to prune:\n\n%s" % tips_to_prune
    
    #TODO: This should be handled by the exclude_tip function (currently in make_test_trees.py)
    #(it has better error handling)
    for tip_name in tips_to_prune:
        tip = new_tree.getNodeMatchingName(tip_name)
        if tip.Parent is not None:
            removal_ok = tip.Parent.remove(tip)
        else:
            removal_ok = False
        new_tree.prune()
    
    return new_tree

def print_node_summary_table(input_tree):
    """Print a summary of the name,children,length, and parents of each node"""
    for node in input_tree.postorder():
        if node.Parent:
            parent_name = node.Parent.Name
        else:
            parent_name = None
        yield "\t".join(map(str,[node.Name,len(node.Children),node.Length,parent_name]))
 

def add_to_filename(filename,new_suffix,delimiter="_"):
    """Add to a filename, preserving the extension"""
    filename, ext = splitext(filename)
    new_filename = delimiter.join([filename,new_suffix])
    return "".join([new_filename,ext])
 

def make_id_mapping_dicts(tree_to_trait_mappings):
    """Generates trait_to_tree, tree_to_trait mapping dictionaries from a list of mapping tuples

    mappings -- in the format tree_id, trait_id
    
    """
    trait_to_tree_mapping_dict = {}
    tree_to_trait_mapping_dict = {}
    
    for tree_id,trait_id in tree_to_trait_mappings:
        trait_to_tree_mapping_dict[trait_id] = tree_id
        tree_to_trait_mapping_dict[tree_id] = trait_id

    return trait_to_tree_mapping_dict, tree_to_trait_mapping_dict


def remap_trait_table_organisms(trait_table_lines,trait_to_tree_mapping_dict,\
  default_entry="not_mapped",input_delimiter="\t",output_delimiter="\t"):
    """Yield trait table lines with organism ids substituted using the mapping dict
    
    trait_table_lines -- tab-delimited lines of organism\tcounts (each count gets
      a column).   

    default_entry -- the value that trait table ids should be set to when they can't
    be mapped to ids
    """
 
    #TODO: should just pass list of fields to avoid this nonsense
    
    for line in trait_table_lines:
        if line.startwith("#"):
            continue
        fields = line.strip().split(input_delimiter)

        fields[0] = trait_to_tree_mapping_dict.get(fields[0],default_entry) 
        yield output_delimiter.join(fields)+"\n" 
