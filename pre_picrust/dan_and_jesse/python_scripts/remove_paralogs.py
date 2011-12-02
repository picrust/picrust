from cogent.seqsim.tree import RangeNode
from cogent.parse.tree import DndParser
from collections import defaultdict
from random import choice
from sys import argv
import os
from os import listdir,mkdir

"""Attribution note:  I believe this script was written by some combination
of Julia Goodrich and Daniel McDonald for a differnet project.  
Dan Knights and I were using the remove paralogs
feature to clean up trees built with all 16S copies."""



class SimTreeError(Exception):
    pass


class SimNode(RangeNode):
    def __init__(self, Name=None, LeafRange=None, Id=None, Children=None, \
            Parent=None, NameLoaded=True, Params=None):
        """Returns a new SimNode object.

        Name: text label
        LeafRange: range of this node's leaves in array. last = index+1
        Id: index of this node in array.
        Children: list of Node objects that are this node's direct 
                  children.
        Parent: Node object that is this node's parent.
        NameLoaded: From cogent.core.tree.TreeNode, undocumented.
        Params are the parameters for the nde such as length
        """
        super(SimNode, self).__init__(Name=Name, Children=Children,\
            Parent=Parent, NameLoaded=NameLoaded)
        if Params is None:
            Params = {}
        self.params = Params
        self.LeafRange = LeafRange
        self.Id = Id

    
    def _getSubTree(self, included_names, constructor=None):
        """An equivalent node with possibly fewer children, or None
            this is an iterative version of that in cogent.core.tree.py
        """
        nodes_stack = [[self, len(self.Children)]]
        result = [[]]

        # Renumber autonamed edges
        if constructor is None:
            constructor = self._default_tree_constructor()

        while nodes_stack:
            top = nodes_stack[-1]
            top_node, num_unvisited_children = top
            if top_node.Id in included_names:
                result[-1].append(top_node.deepcopy(constructor=constructor))
                nodes_stack.pop()
            else:
                #check the top node, any children left unvisited?
                if num_unvisited_children: #has any child unvisited
                    top[1] -= 1  #decrease the #of children unvisited
                    next_child = top_node.Children[-num_unvisited_children]
                    # - for order
                    #pre-visit
                    nodes_stack.append([next_child, len(next_child.Children)])
                    if len(next_child.Children) > 0:
                        result.append([])
                else:
                    node = nodes_stack.pop()[0]
                    children = result[-1]
                    children =[child for child in children if child is not None]
                    if len(top_node.Children) == 0:
                        new_node = None
                    elif len(children) == 0:
                        result.pop()
                        new_node = None
                    elif len(children) == 1:
                        result.pop()
                        # Merge parameter dictionaries by adding lengths and
                        # making weighted averages of other parameters.  This
                        # should probably be moved out of here into a
                        # ParameterSet class (Model?) or tree subclass.
                        params = {}
                        child = children[0]
                        
                        if node.Length is not None and child.Length is not None:
                            shared_params = [n for (n,v) in node.params.items()
                                if v is not None
                                and child.params.get(n) is not None
                                and n is not "length"]
                            length = node.Length + child.Length
                            if length:
                                params = dict([(n,
                                        (node.params[n]*node.Length +
                                        child.params[n]*child.Length) / length)
                                    for n in shared_params])
                            params['length'] = length
                        new_node = child
                        new_node.params = params
                    else:
                        result.pop()
                        new_node = constructor(node, tuple(children))
                    if len(result)>0:
                        result[-1].append(new_node)
                    else:
                        return new_node

    def getSubTree(self, name_list):
        """A new instance of a sub tree that contains all the otus that are
        listed in name_list.
        just a small change from that in cogent.core.tree.py so that the root
        node keeps its name
        """
        new_tree = self._getSubTree(name_list)
        if new_tree is None:
            raise SimTreeError, "no tree created in make sub tree"
        elif new_tree.istip():
            # don't keep name
            new_tree.params = self.params
            new_tree.Length = self.Length
            return new_tree
        else:
            new_tree.Name = self.Name
            new_tree.NameLoaded = self.NameLoaded
            new_tree.params = self.params
            new_tree.Length = self.Length
            # keep unrooted
            if len(self.Children) > 2:
                new_tree = new_tree.unrooted()
            return new_tree
 

    def propagateAttrUp(self, attr, overwrite=False):
        node = self
        while node.Parent is not None and not hasattr(node.Parent, attr):
            setattr(node.Parent, attr, getattr(node,attr))
            node = node.Parent


if __name__ == "__main__":
    #high_directory = argv[1]
    #directory = os.path.join(high_directory,"genes_by_ko")
    #print directory
    #new_dir = directory.split('/')[:-1]
    #new_dir.append("gene_trees_no_paralogs_smaller")
    #new_dir = '/'.join(new_dir)
    #print new_dir
    #try:
    #    mkdir(new_dir)
    #except OSError:
    #    pass
    #for f in listdir(directory):
    tree = DndParser(open(argv[1]).read(), constructor = SimNode)
    lp = tree.getTipNames()
    lp.sort()
    # print "Before tree num", len(lp)
    #print "Before tree set",len(set(lp))
    tree.assignIds()
    tips = defaultdict(list)
    for n in tree.tips():
        new_name = n.Name.split('_')[0].strip("'")
        n.Name = new_name
        tips[n.Name].append(n)
    to_keep = []
    for name, n in tips.items():
        #to_keep.append(choice(n).Id)
        if name.strip("'").startswith("METSMI"):
            if len(n) > 1:
                print "we have a problem"
            to_keep.append(n[0].Id)
            continue
        n_by_length = []
        n_by_length_dict = defaultdict(list)
        for a in n:
            n_by_length.append((a.Length,a))
            n_by_length_dict[a.Length].append(a)
        #print n_by_length
        min_len = min(n_by_length)[0]
        tips_for_min = n_by_length_dict[min_len]
        #print tips_for_min
        to_keep.append(choice(tips_for_min).Id)
    new_tree = tree.getSubTree(to_keep)
    new_f = open(argv[1] + '.noparalog','w')
    new_f.write(new_tree.getNewick(with_distances = True))
    new_f.close()

