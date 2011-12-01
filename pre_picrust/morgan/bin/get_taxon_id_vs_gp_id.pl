#!/usr/bin/perl

use warnings;
use strict;
use autodie;

use lib "/home/mlangill/Dropbox/projects/";

#we need to use the Search library (this also imports GenomeProject,Replicon, and Gene libs)
use MicrobeDB::Search;

my $version_id=25;

#load in the pfam predictions dataset so we don't return 16s sequences that we don't have data for
my $data_file="../../data/pfam_vs_gp_id.txt";
open(my $DATA,'<',$data_file);
my $header=<$DATA>;
chomp($header);
my ($blah,@gp_ids)=split(/\s+/,$header);

my %data_gp_ids;
foreach(@gp_ids){
    $data_gp_ids{$_}=1;
}

#Create the search object.
my $so = new MicrobeDB::Search();

my @gpos = $so->object_search(new MicrobeDB::GenomeProject(version_id=>$version_id));

foreach my $gpo (@gpos) {
    
    my $gpv_id=$gpo->gpv_id();
    my $gp_id=$gpo->gp_id();
    
    #don't bother with genomes where we don't have pfam predictions
    unless($data_gp_ids{$gp_id}){
	#warn "No pfam predictions for genome project $gpv_id";
	next;
    }
    my $taxon_id= $gpo->taxon_id();

    if(defined($taxon_id)){

	print $gp_id,"\t",$taxon_id,"\n";
    }else{
	warn "No taxon_id in MicrobeDB for gp_id: $gp_id";
    }
}
