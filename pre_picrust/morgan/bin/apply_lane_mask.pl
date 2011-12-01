#!/usr/bin/env perl

use warnings;
use strict;
use autodie;
use Bio::AlignIO;
use File::Basename;
use Cwd 'abs_path';

my $abs_dir=dirname(abs_path($0)).'/';
$abs_dir=dirname($abs_dir).'/';

my $in = Bio::AlignIO->newFh('-format' =>'Fasta','-fh'=>\*ARGV);
my $out = Bio::AlignIO->newFh('-format'=>'Fasta','-fh'=>\*STDOUT);

my $lanemask_file =$abs_dir. 'data/lanemask_in_1s_and_0s';

open (my $LANE,'<',$lanemask_file);

my @lines=<$LANE>;
chomp(@lines);
my $lane=join('',@lines);
my @lanemask =split('',$lane);

my $start;
my @bad_columns;
#create an array of array refs that signify start and stop positions of columns to take out
for my $i (0..$#lanemask){
    if(!defined($start) && ! $lanemask[$i]){
	#found the start of a new area to remove
	$start=$i;
    }elsif(defined($start) && $lanemask[$i]){
	#within a remove area and found end of remove area
	push @bad_columns,[$start,$i-1];
	undef $start;
    }
}
#have to handle special case of 0's until end
if(defined($start)){
    push @bad_columns,[$start,$#lanemask];
}

my $align = <$in>;
my $new_align = $align->remove_columns(@bad_columns);
print $out $new_align;



