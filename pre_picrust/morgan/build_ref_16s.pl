#!/usr/bin/env perl
#perldoc build_ref_16s.pl

use warnings;
use strict;
use File::Basename;
use Log::Log4perl;
use Getopt::Long;
use Pod::Usage;
use Cwd 'abs_path';

#Set up options
my %opt=();
GetOptions (\%opt,'Force','help') or pod2usage(2);
pod2usage(-verbose=>2) if exists $opt{'help'};

my $fasta_file =$ARGV[0];

pod2usage($0.': You must provide a fasta file as input') if !defined($fasta_file)|| ! -e $fasta_file;

my $abs_dir=dirname(abs_path($0)).'/';

# Set up the logger             
my $logger_cfg = $abs_dir."logger.conf";
Log::Log4perl::init($logger_cfg);
my $logger = Log::Log4perl->get_logger;

my $name=fileparse($fasta_file);
my $tmp_dir = $abs_dir.'ref_trees/'.$name.'/';
$logger->info("Going to store all files in $tmp_dir");
system("mkdir -p $tmp_dir");
if ($? == -1) {
	$logger->fatal("failed to execute: $!");
	exit;
}


my $template_file= $abs_dir."data/core_set_aligned.fasta.imputed";


#create alignment using pynast
my $pynast_alignment_file=$tmp_dir .'pynast_alignment.fa';
my $pynast_fail_file=$tmp_dir.'pynast_fail.fa';
my $pynast_log_file=$tmp_dir.'pynast_log.txt';
my $pynast_cmd="pynast -g $pynast_log_file -f $pynast_fail_file -t $template_file -a $pynast_alignment_file -i $fasta_file";
$logger->info("Creating alignment of reference sequences using pynast with green genes template");
$logger->debug($pynast_cmd);
system($pynast_cmd);
if ($? != 0) {
	$logger->fatal("failed to execute: $!");
	exit;
}

#trim alignment using lane mask
my $pynast_trimmed_alignment = $tmp_dir.'pynast_trimmed_alignment.fa';
my $trim_cmd = $abs_dir."bin/apply_lane_mask.pl $pynast_alignment_file | ".$abs_dir."bin/fix_bioperl_msa_ids.pl > $pynast_trimmed_alignment";
$logger->info("Trimming alignment using lane mask");
$logger->debug($trim_cmd);
system($trim_cmd);
if ($? != 0) {
	$logger->fatal("failed to execute: $!");
	exit;
}


#convert alignment from fasta to phylip
my $pynast_trimmed_alignment_phylip = $tmp_dir.'pynast_trimmed_alignment.phylip';
my $convert_cmd = $abs_dir."bin/convert_alignment.pl $pynast_trimmed_alignment >$pynast_trimmed_alignment_phylip";
$logger->info("Converting alignment to phylip format for use with raxml");
$logger->debug($convert_cmd);
system($convert_cmd);
if ($? != 0) {
	$logger->fatal("failed to execute: $!");
	exit;
}


#create tree using raxml
my $raxml_output_name="16s";
my $raxml_tree_file = $tmp_dir .'RAxML_result.'.$raxml_output_name;
my $raxml_stats_file = $tmp_dir.'RAxML_info.'.$raxml_output_name;
if(-e $raxml_tree_file && ! $opt{'Force'}){
    $logger->info("Skipping creation of RAxML tree since tree file already exists");
}else{
    my $raxml_cmd="raxmlHPC -T 2 -m GTRGAMMA -w $tmp_dir -n $raxml_output_name -s $pynast_trimmed_alignment_phylip";
    $logger->info("Creating tree using RAxML");
    $logger->debug($raxml_cmd);
    system($raxml_cmd);
    if ($? != 0) {
	$logger->fatal("failed to execute: $!");
	exit;
    }
}

__END__

=head1 Name

build_ref_16s.pl - Creates a reference 16S phylogenetic tree from a set of 16S sequences 

=head1 USAGE

build_ref_16s.pl [OPTIONS] <16S_FASTA_FILE>

E.g.:

build_ref_16s.pl 16s_seqs_in_seed.fa

=head1 OPTIONS

=over 4

=item B<-F, --Force>

Forces all steps of tree constuction even if output files already exist. 

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<build_ref_16S.pl> starts by aligning the 16S sequences to the green genes reference alignment using pynast. Then the lane mask is applied to the alignment to remove "bad" columns. The alignment is then converted from FASTA format to phylip format so that it can be used as input for RAxML tree construction.

=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

=head1 DATE

08-Dec-2011

=cut

