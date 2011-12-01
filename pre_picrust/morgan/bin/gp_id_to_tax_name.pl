#!/usr/bin/env perl

use warnings;
use strict;

use lib "/home/mlangill/Dropbox/projects/";
use MicrobeDB::Search;

my $so= new MicrobeDB::Search();

my @gpos= $so->object_search(new MicrobeDB::GenomeProject(version_id=>25));

foreach my $gpo (@gpos){
    my $gp_id=$gpo->gp_id();
    my $org_name=$gpo->gp_id();
    print $gp_id,"\t",$org_name,"\n";
}
