#!/usr/bin/perl

use warnings;
use strict;

my $pfam_meta_file='data/pfam_metadata.txt';
my $abundance_min_cutoff='3';
my $abundance_max_cutoff='5000';

open(my $META,'<',$pfam_meta_file) || die "Can't open $pfam_meta_file for reading: $!";

#discard header
<$META>;

my %abundance;
while(<$META>){
    chomp;
    my ($pfam_id,$abundance_count)=(split);
    $abundance{$pfam_id}=$abundance_count;
}

my $header=<>;
print $header;
while(<>){
    my ($pfam_id)=split;
    die "Can't find pfam_id: $pfam_id in $pfam_meta_file" unless exists($abundance{$pfam_id});
    if($abundance{$pfam_id}>=$abundance_min_cutoff && $abundance{$pfam_id}<=$abundance_max_cutoff){
	print $_;
    }
}
