#!/usr/bin/env perl

#pulls out the ids from a fasta file and outputs them to stdout

use warnings;
use strict;
use Bio::SeqIO;

my $IN = Bio::SeqIO->newFh('-format' =>'Fasta','-fh'=>\*ARGV);

while(<$IN>){
    my $id = $_->display_id;
    my $seq = new Bio::Seq('-display_id'=>$id,'-seq'=>$_->seq());
    print $id,"\n";
}
