#!/usr/bin/perl

#Copyright (C) 2011 Morgan G.I. Langille
#Author contact: morgan.g.i.langille@gmail.com

#This file is part of MicrobeDB.

#MicrobeDB is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#MicrobeDB is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with MicrobeDB.  If not, see <http://www.gnu.org/licenses/>.

#Example of how to use the search api to get information from microbedb using an object as the search field
#Searchable objects are:
#GenomeProject, Replicon, Gene, Version, or UpdateLog

#See table_search_example.pl, if you want to do a simple search on a mysql db table that is not part of the microbedb api  

#This script retrieves all annotated 16s genes and outputs them in fasta file format.

use warnings;
use strict;
use autodie;

use lib "/home/mlangill/Dropbox/projects/";

#we need to use the Search library (this also imports GenomeProject,Replicon, and Gene libs)
use MicrobeDB::Search;

my $version_id=25;

#load in the pfam predictions dataset so we don't return 16s sequences that we don't have data for
my $data_file="../data/pfam_vs_gp_id.txt";
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
    
    unless($data_gp_ids{$gp_id}){
	warn "No pfam predictions for genome project $gpv_id";
	next;
    }

    my @genes = $so->object_search(new MicrobeDB::Gene(gpv_id=>$gpv_id,gene_type=>'rRNA'));
    unless (@genes){
	warn "No rRNA genes for genome project $gpv_id";
	next;
    }
    
    my $gene = get_single_16s_gene(\@genes);
    unless(defined($gene)){
	warn "No 16s gene for genome project $gpv_id";
	next;
    }

    unless($gene->gene_length() > 1000){
	warn "16s gene is too short for genome project $gpv_id";
	next;
    }
    my $seq = $gene->gene_seq();

    #reverse complement the gene if on negative strand
    if($gene->gene_strand()==-1){
	$seq= revdnacomp($seq);
    }
   
    #print out the gene in fasta format
    print ">$gp_id\n";
    print $seq . "\n";
}

sub revdnacomp {
  my $dna = shift;
  my $revcomp = reverse($dna);

  $revcomp =~ tr/ACGTacgt/TGCAtgca/;

  return $revcomp;
}

sub get_single_16s_gene{
    my $rRNAs=shift;
    my @genes;
    foreach my $rRNA (@$rRNAs){
	my $gene_product = $rRNA->gene_product();
	if(defined($gene_product) && $gene_product =~ /16s/i){
	    push @genes,$rRNA;
	}
    }
    #return nothing if no 16s genes found
    return unless @genes;
    #return single 16s gene if only 1 found
    return $genes[0] if scalar(@genes)==1;

    #if 16s gene names are used then pick the one with name rrsa (case insensitive)
    foreach my $gene(@genes){
	my $name=$gene->gene_name();
	if( defined($name) && $name =~/rrsa/i){
	    return $gene;
	}
    }
    #pick the longest 16s gene (if same length first one is chosen
    my $longest=0;
    my $longest_16s;
    foreach my $gene(@genes){
	if($gene->gene_length() > $longest){
	    $longest_16s=$gene;
	}
    }
    return $longest_16s;
}
