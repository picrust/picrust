#!/usr/bin/env perl
#converts a multiple alignment file in FASTA format to phylip format using bioperl
use warnings;
use strict;

use Bio::AlignIO;

my $in = Bio::AlignIO->newFh('-format' =>'Fasta','-fh'=>\*ARGV);
my $out = Bio::AlignIO->newFh('-format'=>'Phylip','-fh'=>\*STDOUT);
my $aln=<$in>;
print $out $aln;



