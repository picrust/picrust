#!/usr/bin/env perl

use warnings;
use strict;
use autodie;
use Getopt::Long qw(:config no_ignore_case);
use Log::Log4perl;
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;

my $abs_dir=dirname(abs_path($0)).'/';
$abs_dir=dirname($abs_dir).'/';

#Set up options
my %opt=();
GetOptions (\%opt,'input=s','help') or pod2usage(2);
pod2usage(-verbose=>2) if exists $opt{'help'};

my $ref_name=$opt{'input'};

my $ref= $abs_dir.$ref_name;

##build reference tree
my $build_ref_cmd=$abs_dir."build_ref_16s.pl $ref";
#system($build_ref_cmd);


#get ids from seq file
my $grep_cmd='grep ">" '.$ref.' | cut -c 2-';
my @ids=`$grep_cmd`;
chomp(@ids);


my $orig_ref_dir=$abs_dir.'ref_trees/'.$ref_name;
my $orig_ref_tree= $orig_ref_dir.'/RAxML_result.16s';
my $orig_ref_align=$orig_ref_dir.'/pynast_trimmed_alignment.fa';

#for each genome
my $count=0;
foreach my $id (@ids){
    $count++;
    print "$count\n";
    my $new_ref_name=$ref_name.'_minus_'.$id;
    my $new_ref_dir=$abs_dir.'alignments/';
    system("mkdir -p $new_ref_dir");
    my $new_ref_align=$new_ref_dir.$new_ref_name;
    my $new_phylip_align=$new_ref_align.".phylip";

    next if -e $new_phylip_align;

    #remove the 16s from the alignment
    my $remove_align_cmd=$abs_dir."bin/remove_seq_from_alignment.pl $orig_ref_align $id $new_ref_align";
    system($remove_align_cmd);
  
    #convert it to phylip
    my $convert_cmd = $abs_dir."bin/convert_alignment.pl $new_ref_align >$new_phylip_align";
    system($convert_cmd);

    system("rm $new_ref_align");

}
