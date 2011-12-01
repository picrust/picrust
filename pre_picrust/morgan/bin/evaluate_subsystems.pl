#!/usr/bin/perl

use warnings;
use strict;
use autodie;
use List::Util qw(sum);
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

my $abs_dir=dirname(abs_path($0)).'/';
$abs_dir=dirname($abs_dir).'/';

my $method;
my $func;
GetOptions("method=s" =>\$method, "func=s"=>\$func);

#load in the subsystems with hiearchael information
my $subs_with_cat_file=$abs_dir.'data/subsystem_classes.txt';
open(my $SUB_INFO,'<',$subs_with_cat_file);

my %sub_info;
my %sub_cat;
while(<$SUB_INFO>){
    chomp;
    my($id,$desc,$cat,$subcat)=split(/\t/,$_);
    $sub_info{$id}=$cat;
    $sub_cat{$id}=$cat.'_'.$subcat;
}

#load in the accuracy information
my $sub_accuracy_file=$abs_dir.'results/subsystem_neighbour_precision_func_accuracy.txt';
open(my $ACCURACY,'<',$sub_accuracy_file);
#throw away header
<$ACCURACY>;
my %accuracy;
while(<$ACCURACY>){
    chomp;
    my ($id,$acc)=split;
    $accuracy{$id}=$acc;
}

my %category_acc;
foreach my $id(keys(%sub_info)){
    push(@{$category_acc{$sub_info{$id}}},$accuracy{$id});
}

my %category_avg;
my $delim="\t";
#print header
print join($delim,'category','num_of_ss_within_cat','average_accuracy'),"\n";
foreach my $cat (keys(%category_acc)){
    print join($delim, $cat,scalar(@{$category_acc{$cat}}),sum(@{$category_acc{$cat}})/scalar(@{$category_acc{$cat}})),"\n";
}
