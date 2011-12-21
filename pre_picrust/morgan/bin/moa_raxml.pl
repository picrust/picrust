#!/usr/bin/env perl
#This runs RaXML on all alignments in the input directory.

use warnings;
use strict;

use File::Basename;
use File::Copy;
use Time::HiRes qw(usleep);
use POSIX qw/ceil/;

use Cwd;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;

my $abs_dir=dirname(abs_path($0)).'/';
$abs_dir=dirname($abs_dir).'/';

my $filter='*';
my $in_dir;
my $out_dir;
my $database;
my $pipe_opt;
my $gzip;
my $help;
GetOptions(
    "in_dir=s"=>\$in_dir,
    "out_dir=s"=>\$out_dir,
    "filter=s" => \$filter,
    "database=s" =>\$database,
    "pipe_opt=s" =>\$pipe_opt,
    "zip" =>\$gzip,
    "help"           => \$help
    );

my $usage = 'moa_raxml.pl [-h -z -p <input pipe option> -f <file input filter> ] -i <input_directory> -o <output_directory>' . "\n\n"; 

if($help || !defined($in_dir) || !defined($out_dir)){
    print $usage;
    exit;
}

unless(-d $in_dir){
    die "$in_dir is not a directory";
}
my @input_files=glob($in_dir.'/'.$filter);


unless(@input_files){
    die "No input files were found!";
}

my $num_of_nodes=300;
my $num_files_limit=ceil(scalar(@input_files)/$num_of_nodes);


my $name='tree';
my $node_scratch = '/scratch/mlangill';
my $qjob_scratch = '/home/mlangill/tmp/scratch/';

##  All your wrappers are written in directory scratch
## Jobs are also run from the same directory [Remember -cwd option]
## Create a directory scracth if not already created

my $command = "rm -rf $qjob_scratch" .'*';
system($command);

system("mkdir -p $out_dir");

my $cdir = getcwd();
#print "$cdir\n";

my $dir_count=0;
my $file_count=0;
my $in_name="moa_".$dir_count;
my $tmp_dir=$in_dir.'/'.$in_name;
system("mkdir -p $tmp_dir");

foreach my $input_file (@input_files){

    $file_count++;    
    #move the file to tmp dir
    move($input_file,$tmp_dir) or die "Move failed: $!";
    
    if($file_count > $num_files_limit){
	$file_count=0;
	$dir_count++;
	## Create Wrapper File
    
	my $qjob_file = qjob($tmp_dir,$out_dir, $in_name);
    
	## Go to the scratch directory
	chdir $qjob_scratch;
    
	### Submit job (-V imports environment variables)
	#$command = "qsub -V $qjob_file";
	$command = "qsub $qjob_file";
	system($command);

	chdir($cdir);
	## Give the poor computer some rest
	usleep(1000*500);

	#create a new directory
	$in_name="moa_".$dir_count;
	$tmp_dir=$in_dir.'/'.$in_name;
	system("mkdir -p $tmp_dir");
    }

}

#launch the last job

my $qjob_file = qjob($tmp_dir,$out_dir, $in_name);

## Go to the scratch directory
chdir $qjob_scratch;

### Submit job (-V imports environment variables)
#$command = "qsub -V $qjob_file";
$command = "qsub $qjob_file";
system($command);




sub qjob{
    my ($input_dir,$out_dir,$input_name) = @_;
        
    my $qjob_file = $qjob_scratch . "ML_".$input_name.".sh";
    print "$qjob_file\n";
    open(FOUT,'>',$qjob_file) || die "Can't open: $qjob_file for writing. $!";
    print FOUT "#!/bin/bash\n";
    print FOUT "##\n";
    print FOUT "#\$ -cwd\n";
    #moa requires a wall time, lets use 96 hours
    print FOUT "#\$ -l h_rt=96:00:00\n";
    ###
    print FOUT "#\$ -S /bin/bash\n";
    print FOUT "#\$ -o $out_dir/$input_name.out\n";
    print FOUT "#\$ -e $out_dir/$input_name.err\n";
    print FOUT "hostname\n";
   
    my $command .= $abs_dir."bin/run_raxml.pl -i $input_dir -o $out_dir";
    #command to be run
    print FOUT "$command\n";
    
    close(FOUT);
    return $qjob_file;
}
