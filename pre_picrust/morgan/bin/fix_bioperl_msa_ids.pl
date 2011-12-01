#!/usr/bin/env perl

use warnings;
use strict;
use autodie;
use Bio::SeqIO;

my $IN = Bio::SeqIO->newFh('-format' =>'Fasta','-fh'=>\*STDIN);
my $OUT = Bio::SeqIO->newFh('-format'=>'Fasta','-fh'=>\*STDOUT);

while(<$IN>){
    my $id = $_->display_id;
    ($id)=split(/\//,$id);
    my $seq = new Bio::Seq('-display_id'=>$id,'-seq'=>$_->seq());
    print $OUT $seq;
}



