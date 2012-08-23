#!/usr/bin/env perl
#converts a multiple alignment file in FASTA format to phylip format using bioperl
use warnings;
use strict;

use Bio::AlignIO;

my $in = Bio::AlignIO->newFh('-format' =>'Fasta','-fh'=>\*ARGV);
my $aln=<$in>;
my $len_longest_id=$aln->maxdisplayname_length();
my $out = Bio::AlignIO->newFh('-format'=>'Phylip','-fh'=>\*STDOUT, '-idlength'=>$len_longest_id);
print $out $aln;



