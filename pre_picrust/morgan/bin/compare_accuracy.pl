#!/usr/bin/env perl

use warnings;
use strict;
use autodie;
use List::Util qw(sum);
use Getopt::Long;
use File::Basename;

my $first_file= $ARGV[0];
my $second_file=$ARGV[1];

my $one_name=fileparse($first_file);
my $two_name=fileparse($second_file);

open(my $ONE,'<',$first_file);
open(my $TWO,'<',$second_file);

my %one;
while(<$ONE>){
    chomp;
    my ($id,$value)=split;
    $one{$id}=$value;
}

my %two;
while(<$TWO>){
    chomp;
    my ($id,$value)=split;
    $two{$id}=$value;
}

my $same=0;
my $better_one=0;
my $better_two=0;
my $better_one_sum=0;
my $better_two_sum=0;

foreach my $id(keys %one){
    unless(exists($two{$id})){
	warn "ID: $id not found in file: $two_name";
	next;
    }
    my $diff_value= $one{$id}-$two{$id};
    if($diff_value>0){
	$better_one++;
	$better_one_sum+=$diff_value;
    }elsif($diff_value<0){
	$better_two++;
	$better_two_sum+=$diff_value;
    }else{
	$same++;
    }
    
    print $id,"\t",$diff_value,"\n";

}

my $better_one_avg = $better_one_sum/$better_one;
my $better_two_avg = $better_two_sum/$better_two;

print STDERR "----SUMMARY----\n";
print STDERR "Same count: $same\n";
print STDERR "Better in $one_name count: $better_one\n";
print STDERR "Better in $one_name average: $better_one_avg\n";
print STDERR "Better in $two_name count: $better_two\n";
print STDERR "Better in $two_name average: $better_two_avg\n";
