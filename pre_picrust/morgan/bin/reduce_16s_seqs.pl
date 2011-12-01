#!/usr/bin/perl

use warnings;
use strict;

use Bio::SeqIO;

my $gp_list="../data/gp_list_in_seed.txt";
open(my $LIST,'<',$gp_list);

my %list;
while(<$LIST>){
    chomp;
    $list{$_}=1;
}

my $IN = Bio::SeqIO->newFh(-format=>'Fasta', -fh =>\*ARGV);
my $OUT = Bio::SeqIO->newFh(-format=>'Fasta', -file =>">16s_seqs_in_seed.fa");

while(<$IN>){
    if($list{$_->id()}){
	print $OUT $_;
    }
}
