#!/usr/bin/env perl

use warnings;
use strict;
use File::Basename;
use Cwd 'abs_path';

my $abs_dir=dirname(abs_path($0)).'/';
$abs_dir=dirname($abs_dir).'/';


my $results_dir=$abs_dir.'results/';

system ("mkdir -p $results_dir");

foreach my $func ('pfam','subsystem','role','EC'){
    foreach my $method ('neighbour','random','pic'){   
	foreach my $acc_type ('precision','recall'){
	    my $acc_file = $func.'_'.$method.'_'.$acc_type.'_accuracy.txt';
	    my $cat_cmd = "cat $abs_dir"."tmp/*/$acc_file > ".$results_dir.$acc_file;
	    system($cat_cmd);
	}
    }
}
