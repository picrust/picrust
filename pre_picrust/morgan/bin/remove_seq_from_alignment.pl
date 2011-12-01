#!/usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO;
my $orig_align_file= $ARGV[0];
my $id_to_remove= $ARGV[1];
my $new_align_file=$ARGV[2];


my $IN = Bio::SeqIO->newFh(-file=>$orig_align_file, -format=>'fasta');
my $OUT = Bio::SeqIO->newFh(-file=>">$new_align_file",-format=>'fasta');

while(<$IN>){
    print $OUT $_ unless ($_->id() eq $id_to_remove);
}
