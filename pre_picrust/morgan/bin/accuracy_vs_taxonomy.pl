#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use lib "/home/mlangill/Dropbox/projects/";
use MicrobeDB::Search;
my $taxon_level="phylum";
my $suffix="";
GetOptions("taxon=s"=>\$taxon_level,"suffix=s"=>\$suffix);

my $so= new MicrobeDB::Search();

my @gpos= $so->object_search(new MicrobeDB::GenomeProject(version_id=>25));

my %gp_table;
while(<>){
    chomp;
    my ($gp_id,$accuracy)=split;
    $gp_table{$gp_id}{acc}=$accuracy;
}

my %taxon_count;
foreach my $gp (@gpos){
    my $gp_id=$gp->gp_id();
    my $taxon_high= $gp->phylum();
    my $taxon = $gp->$taxon_level();
    if(!defined($taxon_high) || $taxon_high eq ''){
	$taxon_high ='Other';
    }
    if(!defined($taxon) || $taxon eq ''){
	$taxon ='Other';
    }
    if(exists($gp_table{$gp_id})){
	$gp_table{$gp_id}{high}=$taxon_high;
	$gp_table{$gp_id}{taxon}=$taxon;
	$taxon_count{$taxon_high}{$taxon}++;
    }
}
foreach my $gp_id (keys %gp_table){
    my $taxon = $gp_table{$gp_id}{taxon};
    my $taxon_high=$gp_table{$gp_id}{high};
    my $high_char=substr($taxon_high,0,1);
    my $count = $taxon_count{$taxon_high}{$taxon};
    my $label= join('_',$high_char,$taxon,$suffix,$count);
    
    #print join("\t",$taxon_high,$label,$gp_table{$gp_id}{acc}),"\n";
    
    print join("\t",$count,$gp_table{$gp_id}{acc}),"\n";
    
}

