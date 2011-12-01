#!/usr/bin/perl

use warnings;
use strict;
use autodie;
use File::Basename;
use Cwd 'abs_path';

my $abs_dir=dirname(abs_path($0)).'/';
$abs_dir=dirname($abs_dir).'/';
open my $FH, '<',$abs_dir.'data/gp_id_rep_accnum_lookup.tab';

my %lookup;
while(<$FH>){
    chomp;
    my ($gp_id,$rep_accnum)=split;
    $lookup{$rep_accnum}=$gp_id;
}

while(<>){
    chomp;
    my ($read_id,$ref,$dist)=split(/ /,$_);
    my $gp_id;
    if(exists($lookup{$ref})){
	$gp_id = $lookup{$ref};
    }else{
	warn "Couldn't convert $ref";
	next;
    }
    print join("\t",$read_id,$lookup{$ref},$dist),"\n";
}
