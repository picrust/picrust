#!/usr/bin/env perl

use warnings;
use strict;
use autodie;
use Parallel::ForkManager;
use Log::Log4perl;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

my($help);
my $threads=0;
my $ref_name="16s_seqs_in_seed.fa";
my $figures_flag;
GetOptions ("threads=i" => \$threads, "ref_fasta=s"=>\$ref_name, 'figures!'=>\$figures_flag,"help"=>\$help);

my $usage ="$0 [--threads <# num of proc to use> --ref_fasta <FASTA file with all 16S> --figures --nofigures] \n";
if($help){
    print $usage;
    exit;
}

my $abs_dir=dirname(abs_path($0)).'/';

# Set up the logger             
my $logger_cfg = $abs_dir."logger.conf";
Log::Log4perl::init($logger_cfg);
my $logger = Log::Log4perl->get_logger;

#Allows multiple threads to be used
my $pm = new Parallel::ForkManager($threads);


my $ref= $abs_dir.$ref_name;

##build reference tree
my $orig_ref_dir=$abs_dir.'ref_trees/'.$ref_name;
my $orig_ref_tree= $orig_ref_dir.'/RAxML_result.16s';
my $orig_ref_align=$orig_ref_dir.'/pynast_trimmed_alignment.fa';

my $build_ref_cmd=$abs_dir."build_ref_16s.pl $ref";
system($build_ref_cmd) unless -e $orig_ref_tree;

###!!!create leave one out trees (run commands in "create_leave_one_trees" dir)
#This creates leave one out trees on moa cluster, then copies them into proper ref_trees directories

#get ids from seq file
my $grep_cmd='grep ">" '.$ref.' | cut -c 2-';
my @ids=`$grep_cmd`;
chomp(@ids);

   
my $query_dir = $abs_dir."queries/";
system ("mkdir -p $query_dir");

my @func_types=("subsystem","EC","role","pfam");

#do entire pipeline unless --figures flag is set
unless($figures_flag==1){

#for each genome
foreach my $id (@ids){
    
    my $new_ref_name=$ref_name.'_minus_'.$id;
    my $new_ref_dir=$abs_dir.'ref_trees/'.$new_ref_name;
    my $new_ref_tree= $new_ref_dir.'/RAxML_result.16s';

    #don't continue unless tree is already made (from running seperate scripts as described above)
    next unless -e $new_ref_tree;
    #create the new ref dir
    #system("mkdir -p $new_ref_dir");

    #fork (if --threads option is used)
    my $pid= $pm->start and next;

    $logger->debug("Running leave one out for genome: $id");

    #remove the 16s from the alignment
    my $new_ref_align=$new_ref_dir.'/pynast_trimmed_alignment.fa';
    my $remove_align_cmd=$abs_dir."bin/remove_seq_from_alignment.pl $orig_ref_align $id $new_ref_align";
    system($remove_align_cmd) unless -e $new_ref_align;

    #create a query file with just the single 16s 
    my $query_name=$id.'_query.fa';
    my $query= $query_dir. $query_name;
    my $create_query_cmd=$abs_dir."bin/create_query_seq.pl $ref $id $query";
    system($create_query_cmd) unless -e $query;

    foreach my $func (@func_types){ 
	#Do ancestral state reconstruction
	my $data_file=$abs_dir."data/".$func."_vs_gp_id.txt";
	my $ace_file=$new_ref_dir."/".$func."_ace_counts.txt.gz";
	my $ace_ref_cmd=$abs_dir."bin/find_ace.R $ace_file $data_file $new_ref_tree pic";
	system($ace_ref_cmd) unless -e $ace_file;
    }

    foreach my $func (@func_types){ 
	#Do ancestral state reconstruction using ML
	my $data_file=$abs_dir."data/".$func."_vs_gp_id.txt";
	my $ace_file=$new_ref_dir."/".$func."_REML_counts.txt.gz";
	my $ace_ref_cmd=$abs_dir."bin/find_ace.R $ace_file $data_file $new_ref_tree ML";
	system($ace_ref_cmd) unless -e $ace_file;
	exit;
    }

    

    #Sometimes query sequence is not seen as a valid 16S by pynast. This creates an empty alignment file.
    #Check for this here and just skip the genome for further analysis if empty.
    my $pynast_query_alignment_file = $abs_dir."tmp/".$query_name."/pynast_alignment.fa";
    next if -z $pynast_query_alignment_file;

    #place the query reads using pplacer
    my $pplacer_file=$abs_dir."tmp/".$query_name."/pplacer.tog.tre";
    my $place_reads_cmd=$abs_dir."place_reads.pl -r $new_ref_name $query";       
    system($place_reads_cmd) unless -e $pplacer_file;

    #check here too for empty alignment file
    next if -z $pynast_query_alignment_file;

	
    foreach my $func (@func_types){

	#make predictions using ace
	my $ace_predictions_file=$abs_dir."tmp/".$query_name."/$func"."_ace_predictions.txt";
	my $ace_cmd=$abs_dir."bin/ace.pl -m ace -f $func $query_name $new_ref_name";
	system($ace_cmd) unless -e $ace_predictions_file;
	my $ace_accuracy_file=$abs_dir."tmp/".$query_name."/$func"."_ace_precision_accuracy.txt";
	system($abs_dir."bin/test_accuracy.pl -m ace -f $func $query_name") unless -e $ace_accuracy_file;

	#make predictions using ace with ML
	my $REML_predictions_file=$abs_dir."tmp/".$query_name."/$func"."_REML_predictions.txt";
	my $REML_cmd=$abs_dir."bin/ace.pl -m REML -f $func $query_name $new_ref_name";
	#system($REML_cmd) unless -e $REML_predictions_file;
	my $REML_accuracy_file=$abs_dir."tmp/".$query_name."/$func"."_REML_precision_accuracy.txt";
	#system($abs_dir."bin/test_accuracy.pl -m REML -f $func $query_name") unless -e $REML_accuracy_file;
	

	#make predictions using nearest neighbour

	my $nn_predictions_file=$abs_dir."tmp/".$query_name."/$func"."_neighbour_predictions.txt";
	my $neighbour_cmd=$abs_dir."bin/neighbour.pl -f $func $query_name";
	system($neighbour_cmd) unless -e $nn_predictions_file;
	my $nn_accuracy_file=$abs_dir."tmp/".$query_name."/$func"."_neighbour_precision_accuracy.txt";
	system($abs_dir."bin/test_accuracy.pl -m neighbour -f $func $query_name") unless -e $nn_accuracy_file;
	
	#make predictions using random
	my $random_predictions_file=$abs_dir."tmp/".$query_name."/$func"."_random_predictions.txt";
	my $random_cmd=$abs_dir."bin/random.pl -f $func $query_name";
	system($random_cmd) unless -e $random_predictions_file;
	my $random_accuracy_file=$abs_dir."tmp/".$query_name."/$func"."_random_precision_accuracy.txt";
	system($abs_dir."bin/test_accuracy.pl -m random -f $func $query_name") unless -e $random_accuracy_file;
    }
    $pm->finish;
}
$pm->wait_all_children;

}

#create figures unless the --nofigures option is used.
unless($figures_flag==0){
my $figures_dir=$abs_dir.'figures';
system("mkdir -p $figures_dir");

#get accuracy for each genome (files stored in "/results")
#this calls "bin/violin_plots2.R after concatenating files (results stored in "figures")
my $create_plot_cmd=$abs_dir."bin/create_accuracy_plots.pl";
system($create_plot_cmd);


#get accuracy for each pfam (files stored in "/results")
foreach my $func (@func_types){
    foreach my $method ("neighbour","random","ace"){
	my $pfam_accuracy_cmd=$abs_dir."bin/test_functional_accuracy.pl -m $method -f $func";
	system($pfam_accuracy_cmd);
    }
}

#create plots about each functional accuracy (results stored in "/figures")
foreach my $func (@func_types){
    my $func_plots_cmd=$abs_dir."bin/func_accuracy_plots.R $func";
    system($func_plots_cmd);
}

}
