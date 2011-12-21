#!/usr/bin/env perl

use warnings;
use strict;
use autodie;
use File::Basename;
use Cwd 'abs_path';
use Getopt::Long;

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

my $ref_tree_dir=$abs_dir.'ref_trees/';

#for each genome
foreach my $id (@ids){
    my $new_ref_name=$ref_name.'_minus_'.$id;

    #create the directory
    my $new_ref_dir=$ref_tree_dir.$new_ref_name;
    system("mkdir -p $new_ref_dir");

    #copy the tree and info file to the new directory
    my $tree=$abs_dir.'trees/RAxML_bestTree.'.$new_ref_name.".phylip.tree";
    my $info=$abs_dir.'trees/RAxML_info.'.$new_ref_name.".phylip.tree";
    my $new_tree=$new_ref_dir."/RAxML_result.16s";
    my $new_info=$new_ref_dir."/RAxML_info.16s";

    next if (-e $new_tree && -e $new_info);
    print $id,"\n";

    if(-e $tree && -e $info){
	if(-s $tree > 2000){
	    system("cp $tree $new_tree");
	    system("cp $info $new_info");
	}else{
	    warn "Tree file is not large enough for: $id";
	}
    }else{
	warn "Tree or info file does not exist for: $id";
    }
    
}
