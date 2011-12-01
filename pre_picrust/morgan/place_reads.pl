#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Log::Log4perl;
use Cwd 'abs_path';

my($ref_tree,$help);
GetOptions ("ref_tree=s" => \$ref_tree,"help"=>\$help);

my $fasta_file =$ARGV[0];

my $abs_dir=dirname(abs_path($0)).'/';

my $name=fileparse($fasta_file);

# Set up the logger             
my $logger_cfg = $abs_dir."logger.conf";
Log::Log4perl::init($logger_cfg);
my $logger = Log::Log4perl->get_logger;


my $ref_dir=$abs_dir."ref_trees/".$ref_tree.'/';
my $raxml_tree_file = $ref_dir.'RAxML_result.16s_with_node_labels';
my $raxml_stats_file= $ref_dir.'RAxML_info.16s';
my $ref_alignment_file = $ref_dir.'pynast_trimmed_alignment.fa';

unless (-e $raxml_tree_file && -e $raxml_stats_file && -e $ref_alignment_file){
    $logger->fatal("Your reference files do not exists at: $ref_dir");
    die;
}

my $tmp_dir = $abs_dir.'tmp/'.$name.'/';
system("mkdir -p $tmp_dir");
if ($? != 0) {
	$logger->fatal("failed to execute: $!");
	die;
}


my $template_file= $abs_dir."data/core_set_aligned.fasta.imputed";


#extract metagenomic read ids (these are needed later with pplacer)
my $fasta_ids=$tmp_dir.'read_ids.txt';
my $extract_id_cmd=$abs_dir."bin/extract_read_ids.pl $fasta_file > $fasta_ids";
$logger->info("Extracting metagenomic read ids from query sequence for use my pplacer");
$logger->debug($extract_id_cmd);
system($extract_id_cmd);
if ($? != 0) {
	$logger->fatal("failed to execute: $!");
	die;
}


#create alignment using pynast
my $pynast_alignment_file=$tmp_dir .'pynast_alignment.fa';
my $pynast_fail_file=$tmp_dir.'pynast_fail.fa';
my $pynast_log_file=$tmp_dir.'pynast_log.txt';
my $pynast_cmd="pynast -g $pynast_log_file -f $pynast_fail_file -t $template_file -a $pynast_alignment_file -i $fasta_file";
$logger->info("Creating alignment of query sequences using pynast with green genes template");
$logger->debug($pynast_cmd);
system($pynast_cmd);
if ($? != 0) {
	$logger->fatal("failed to execute: $!");
	die;
}

#check if pynast didn't find any valid 16S sequences to align
if(-z $pynast_alignment_file || !(-e $pynast_alignment_file)){
    $logger->fatal("Pynast did not find any valid sequences to align");
    die;
}

#trim alignment using lane mask
my $pynast_trimmed_alignment = $tmp_dir.'pynast_trimmed_alignment.fa';
my $trim_cmd = $abs_dir."bin/apply_lane_mask.pl $pynast_alignment_file | ".$abs_dir."bin/fix_bioperl_msa_ids.pl > $pynast_trimmed_alignment";
$logger->info("Trimming alignment using lane mask");
$logger->debug($trim_cmd);
system($trim_cmd);
if ($? != 0) {
	$logger->fatal("failed to execute: $!");
	die;
}

#place reads onto reference tree using pplacer
my $pplacer_cmd = "pplacer -t $raxml_tree_file -r $ref_alignment_file -s $raxml_stats_file --out-dir $tmp_dir $pynast_trimmed_alignment";
$logger->info("Placing reads on reference tree using pplacer");
$logger->debug($pplacer_cmd);
system($pplacer_cmd);
if ($? != 0) {
	$logger->fatal("failed to execute: $!");
	die;
}

#convert pplacer output to single tree with reads placed on it using guppy
my $placed_reads_tree=$tmp_dir."pplacer";
my $pplacer_json_file=$tmp_dir."pynast_trimmed_alignment.jplace";
my $guppy_cmd = "guppy tog -o $placed_reads_tree $pplacer_json_file";
$logger->info("Converting pplacer json output to single tree using guppy");
$logger->debug($guppy_cmd);
system($guppy_cmd);
if ($? != 0) {
	$logger->fatal("failed to execute: $!");
	die;
}
