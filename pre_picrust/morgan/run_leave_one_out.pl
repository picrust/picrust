#!/usr/bin/env perl

use warnings;
use strict;
use autodie;
use Parallel::ForkManager;
use Log::Log4perl;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Basename;
use Cwd 'abs_path';

#Set up options
my %opt=();
GetOptions (\%opt,'threads=i','ref_fasta=s','figures!','Force','help') or pod2usage(2);
pod2usage(-verbose=>2) if exists $opt{'help'};

my $ref_name="16s_seqs_in_seed.fa";
if(exists $opt{'ref_fasta'}) {
    $ref_name=$opt{'ref_fasta'};
}

my $threads=$opt{'threads'}||0;

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

#Build leave one out trees
#Do this seperately in case we need to run on moa (compute cluster)
if(0){
    #Build leave one out alignments
    my $leave_out_aligns_cmd= $abs_dir.'bin/create_alignments.pl -i '.$ref_name;
    system($leave_out_aligns_cmd);

    #Build trees (using a single node with 8 processors)
    my $leave_out_trees_cmd=$abs_dir.'bin/run_raxml.pl -p 8 -i '.$abs_dir.'alignments/ -o '.$abs_dir.'trees/';
    system($leave_out_trees_cmd);

    ##OR

    ##Build trees on sge cluster
    my $leave_out_trees_sge_cmd=$abs_dir.'bin/moa_raxml.pl -i '.$abs_dir.'alignments/ -o '.$abs_dir.'trees/';
    #system($leave_out_trees_sge_cmd);
    
    #count the tree files to make sure they are all there
    #ls -1 trees/RAxML_bestTree.16s_seqs.fa_minus_* | wc -l

    #Move the computed trees to their proper location in "ref_trees" directory
    my $leave_out_move_trees_cmd=$abs_dir.'bin/move_trees.pl -i '.$ref_name;
    system($leave_out_move_trees_cmd);
}

#get ids from seq file
my $grep_cmd='grep ">" '.$ref.' | cut -c 2-';
my @ids=`$grep_cmd`;
chomp(@ids);


my $query_dir = $abs_dir."queries/";
system ("mkdir -p $query_dir");

#my @func_types=("subsystem","EC","role","pfam");
my @func_types=("pfam10");

#do entire pipeline unless --figures flag is set
if(!exists $opt{'figures'} ||$opt{'figures'}==0){
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

	my $asr_method='pic';
	#my $asr_method='REML';
	#my $asr_method='ML';

	#Do ancestral state reconstruction
	foreach my $func (@func_types){ 
	    my $ace_ref_cmd=$abs_dir."build_asr.pl --func $func -m $asr_method -r $new_ref_name";
	    system($ace_ref_cmd);
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
	    my $ace_cmd=$abs_dir."make_predictions.pl -m pic --func $func -q $query_name -r $new_ref_name";
	    system($ace_cmd);

	    my $ace_accuracy_file=$abs_dir."tmp/".$query_name."/$func"."_pic_precision_accuracy.txt";
	    system($abs_dir."bin/test_accuracy.pl -m pic --func $func $query_name") unless -e $ace_accuracy_file;

	    #make predictions using ace with ML
	    my $REML_cmd=$abs_dir."make_predictions.pl -m reml --func $func -q $query_name -r $new_ref_name";
	    #system($REML_cmd) unless -e $REML_predictions_file;
	    my $REML_accuracy_file=$abs_dir."tmp/".$query_name."/$func"."_reml_precision_accuracy.txt";
	    #system($abs_dir."bin/test_accuracy.pl -m REML -f $func $query_name") unless -e $REML_accuracy_file;
	    

	    #make predictions using nearest neighbour
	    my $neighbour_cmd=$abs_dir."make_predictions.pl -m neighbour --func $func -q $query_name -r $new_ref_name";
	    system($neighbour_cmd);
	    my $nn_accuracy_file=$abs_dir."tmp/".$query_name."/$func"."_neighbour_precision_accuracy.txt";
	    system($abs_dir."bin/test_accuracy.pl -m neighbour --func $func $query_name") unless -e $nn_accuracy_file;
	    
	    #make predictions using random
	    my $random_cmd=$abs_dir."make_predictions.pl -m random --func $func -q $query_name -r $new_ref_name";
	    system($random_cmd);
	    my $random_accuracy_file=$abs_dir."tmp/".$query_name."/$func"."_random_precision_accuracy.txt";
	    system($abs_dir."bin/test_accuracy.pl -m random --func $func $query_name") unless -e $random_accuracy_file;
	}
	$pm->finish;
    }
    $pm->wait_all_children;

}

#create figures unless the --nofigures option is used.
if(!exists $opt{'figures'} || $opt{'figures'}==1){
    my $figures_dir=$abs_dir.'figures';
    system("mkdir -p $figures_dir");

    #get accuracy for each genome (files stored in "/results")
    my $create_plot_cmd=$abs_dir."bin/summarize_accuracy_values.pl";
    system($create_plot_cmd);

    #create violin plots
    my $violin_plot_cmd = $abs_dir."bin/violin_plots2.R";
    system($violin_plot_cmd);


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

__END__

=head1 Name

run_leave_one_out.pl - Calculates accuracy of PICRUST leaving leave one out genome validation method.

=head1 USAGE

run_leave_one_out.pl [OPTIONS] -r=<reference 16S sequences>

E.g.:

run_leave_one_out.pl -r=16s_seqs_in_seed.fa


=head1 OPTIONS

=over 4

=item B<-r, --ref_fasta (REQUIRED)>

Specify the fasta file containing 16S sequences for all completed genomes. 

=item B<-t, --threads>

Pipeline can be multi-threaded. Set this option to number of processors to use. (Default is 1).

=item B<-f, --figures>

Skips past all parts of the pipeline except for generating figures.

=item B<--nofigures>

Runs the pipeline, but does not generate figures.

=item B<-F,--Force>

Forces all steps even if output files already exist (off by default) 

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<run_leave_one_out.pl> calculates the accuracy of the PI-CRUST software using various functional databases and difference methods for prediction (ASR, NN, random). The pipeline removes one genome from the reference genome dataset and treats it as a unknown query (metagenomic) 16S sequence. The pipeline makes predictions on this genome and the predicted functional abundances are compared to the known annotated functional abundances for that genome. This is repeated for every reference genome in the dataset (>1000 times). 

After the pipeline is run, several figures are created to describe the accuracy of PI-CRUST with respect to the genome, distance to nearest neighbour (NN), or on the different functional classes.

=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

=head1 DATE

08-Dec-2011

=cut

