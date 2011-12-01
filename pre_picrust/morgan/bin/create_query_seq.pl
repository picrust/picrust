#!/usr/bin/perl

use warnings;
use strict;

use Bio::SeqIO;

my $orig_seq_file=$ARGV[0];
my $id = $ARGV[1];
my $new_seq_file=$ARGV[2];

my $IN= Bio::SeqIO->newFh(-file=>$orig_seq_file,-format=>'fasta');
my $OUT= Bio::SeqIO->newFh(-file=>">$new_seq_file",-format=>'fasta');

while(<$IN>){
    if($_->id() eq $id){
	print $OUT $_;
	last;
    }
}
