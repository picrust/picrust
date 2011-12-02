#!/usr/bin/env python
#file test_exclude_seqs_by_blast.py

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2010, The GenomeDB project"
__credits__ = ["Jesse Zaneveld", "Rob Knight"]
__license__ = "GPL"
_version__ = "1.0-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Pre-release"

"""
Test code for test_parse_db.py.  
"""
from os import remove, system, mkdir
from random import choice
from numpy import array, arange, log, log10
from cogent.util.unit_test import TestCase, main
from character_state_estimation import find_nearest_unambiguous_parent_state, \
  estimate_state_probs, estimate_tip_states, BINARY
from qiime.parse import parse_newick, PhyloNode
from cogent.seqsim.tree import RangeNode  
from cogent.core.alignment import DenseAlignment
from cogent.seqsim.usage import Rates

class EstimateCharacterTests(TestCase):
    def setUp(self):
        """Set up data for tests"""
        # make hypothetical inferred ancestral state alignments
        anc_fna_lines = ['>root', 'TTATATT',\
                            '>E', 'AAATAAT',\
                            '>F', 'TTATTTT' ]

        self.ancAln = DenseAlignment('\n'.join(anc_fna_lines), MolType = BINARY)


        # make tree
        self.balanced_tree = parse_newick("((A:1,B:1)E:10,(C:1,D:1)F:2)root;", RangeNode)
        
        # make tree
        self.short_balanced_tree = \
          parse_newick("((A:0.000001,B:0.000001)E:10,(C:0.000001,D:0.000001)F:2)root;", RangeNode)
 
        test_q = [[ -1.0, 1.0],[1.0,-1.0]]
        self.equal_q = Rates(array(test_q), self.ancAln.Alphabet**2, normalize=False)

    def test_find_nearest_unambiguous_character_state(self):
        """Count bases should count bases in a sequence"""
                
        # make tip alignments
        tip_fna_lines = [   '>A', 'AA?T?A?',\
                            '>B', 'AAATA??',\
                            '>C', 'TTAT?T?',\
                            '>D', 'TTA????' ]
        tip_aln = DenseAlignment('\n'.join(tip_fna_lines))

        # make tree
 
        nearest_states, dist_states = \
            find_nearest_unambiguous_parent_state(self.balanced_tree, self.ancAln, \
                                         tip_aln, ambig_character ="?")
        #print nearest_states
        #print dist_states

    
    def test_estimate_state_probs(self):
        """Estimate_state_probs should predict state given parent probs, distance and q"""
        
        # Equal probabiity 
        dist_to_parent = 0.10
        parent_probs = {"A":0.50,"T":0.50}
        q = self.equal_q
        probs = estimate_state_probs(parent_probs, dist_to_parent, q)
        self.assertFloatEqual(probs["A"], 0.50)
        self.assertFloatEqual(probs["T"], 0.50)


        # Over the long term, parent probs should be
        # effectively irrelevant: only the rates should matter
        dist_to_parent = 100000.0
        parent_probs = {"A":0.01,"T":0.99}
        test_q = [[ -1.0, 1.0],[1.0,-1.0]]
        q = Rates(array(test_q), self.ancAln.Alphabet**2, normalize=False)
        probs = estimate_state_probs(parent_probs, dist_to_parent, q)
        self.assertFloatEqual(probs["A"], 0.50)
        self.assertFloatEqual(probs["T"], 0.50)
    
    def test_estimate_tip_states(self):
        """Estimate tip states should estimate state probabilities given a tree,
        rate_matrix, and ancestral state probabilities"""

        # single position
        anc_state_probs = [{"E":{"A":0.50,"T":0.50}, "F":{"A":0.01,"T":0.99}}]
        tree = self.short_balanced_tree
        res = estimate_tip_states(tree,self.equal_q, anc_state_probs)  
       
        #~ for position in xrange(len(anc_state_probs)):
            #~ for tip in tree.iterTips():
                #~ print "Position:",position, \
                    #~ 'Tip:', tip.Name, 'Length:', tip.Length,':'
                #~ print "ESTIMATE:",res[position][tip.Name]

        self.assertGreaterThan(res[0]['C']['A'], 0.01)
        self.assertLessThan(res[0]['C']['T'], 0.99)
        self.assertFloatEqual(res[0]['A']['A'], 0.5)
        self.assertFloatEqual(res[0]['A']['T'], 0.5)

        # multiple position
        anc_state_probs = \
            [{"E":{"A":0.50,"T":0.50}, "F":{"A":0.01,"T":0.99}}] * 10
        tree = self.short_balanced_tree
        res = estimate_tip_states(tree,self.equal_q, anc_state_probs)  
       
        for position in xrange(len(anc_state_probs)):
            self.assertGreaterThan(res[position]['C']['A'], 0.01)
            self.assertLessThan(res[position]['C']['T'], 0.99)
            self.assertFloatEqual(res[position]['A']['A'], 0.5)
            self.assertFloatEqual(res[position]['A']['T'], 0.5)
 
if __name__=="__main__":
    main()
