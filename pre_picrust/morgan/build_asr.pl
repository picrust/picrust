#!/usr/bin/env perl
#perldoc build_asr.pl

use warnings;
use strict;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use Log::Log4perl;
use Pod::Usage;
use Cwd 'abs_path';

my $abs_dir=dirname(abs_path($0)).'/';

# Set up the logger             
my $logger_cfg = $abs_dir."logger.conf";
Log::Log4perl::init($logger_cfg);
my $logger = Log::Log4perl->get_logger;

#Set up options
my %opt=();
GetOptions (\%opt,'func=s','method=s','ref_tree=s','Force','help') or pod2usage(2);
pod2usage(-verbose=>2) if exists $opt{'help'};

#Mandatory options
foreach my $option ('ref_tree','method','func'){
    pod2usage($0.': You must specify option for --'.$option) unless exists $opt{$option};
}


my $ref_name=$opt{'ref_tree'};
my $ref_dir=$abs_dir.'ref_trees/'.$ref_name;
my $labelled_ref_tree= $ref_dir.'/RAxML_result.16s';

my $ace_file=$ref_dir.'/'.$opt{'func'}.'_'.$opt{'method'}.'_counts.txt.gz';

if(!-e $ace_file || exists $opt{'Force'}){
    $logger->info("Running $opt{'method'} ASR on $ref_name.");
    
    if($opt{'method'} eq 'wagner'){

	my $data_file=$abs_dir."data/".$opt{'func'}."_vs_gp_id.count";
	$logger->logdie("File $data_file does not exist. Did you specify wrong --func ?") unless -e $data_file;
	my $wagner_results = run_count_wagner($data_file,$labelled_ref_tree);
	open(my $OUT_GZ,"|gzip -c > $ace_file");
	print $OUT_GZ join("\n",@$wagner_results),"\n"; 
	close($OUT_GZ);

    }else{
 
	my $data_file=$abs_dir."data/".$opt{'func'}."_vs_gp_id.txt";
	$logger->logdie("File $data_file does not exist. Did you specify wrong --func ?") unless -e $data_file;

	my $ace_ref_cmd=join(" ",$abs_dir.'bin/find_ace.R',$ace_file, $data_file, $labelled_ref_tree, $opt{'method'});
	$logger->debug("Running cmd: $ace_ref_cmd");
	system($ace_ref_cmd);
    }
}else{
    $logger->info("Not running $opt{'method'} ASR on $ref_name becase ASR output file already exists. Use --Force to overwrite existing files.");
}

sub run_count_wagner{
    my($table,$tree)=@_;
    
    my $count_jar=$abs_dir.'bin/Count.jar';
    my $cmd= join(" ",'java -cp',$count_jar,'ca.umontreal.iro.evolution.genecontent.AsymmetricWagner',$tree,$table,'| grep \'# FAMILY\' ');
    my @output=`$cmd`;    
    chomp(@output);

    #process output

    #process header
    my @head_fields=split(/\t/,$output[0]);
    $output[0]=join("\t",@head_fields[2..$#head_fields-4]);

    #process all other fields
    for my $i(1..$#output){
	my @fields=split(/\t/,$output[$i]);
	$output[$i]=join("\t",@fields[1..$#fields-4]);
    }
    return \@output;
}


__END__

=head1 Name

build_asr.pl - Runs user's choice of ancestral state reconstruction method on a given tree 

=head1 USAGE

build_asr.pl [OPTIONS] -f=[pfam|EC|subsystem|role] -m=[pic|ML|REML|wagner] -r=<reference_tree_name>

E.g.:

build_asr.pl -f=pfam -m=pic -r=16s_seqs_in_seed.fa

=head1 OPTIONS

=over 4

=item B<-m, --method>

Defines which asr method to run. Choices include: pic,ML,REML

=item B<-f,--func>

Defines which functional traits are being used for reconstruction. Choices include: subsystem, pfam, EC, roles

=item B<-r, --ref_tree>

Specify reference tree to run ASR on. Note this is just the directory name within (ref_trees).

=item B<-F,--Force>

Forces asr to run even if asr output files already exist(non-default) 

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<build_asr.pl> creates ancestral state reconstruction using one of various methods for a given reference tree.

=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

=head1 DATE

08-Dec-2011

=cut

