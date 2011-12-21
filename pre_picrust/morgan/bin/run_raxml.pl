#!/usr/bin/perl

use warnings;
use strict;
use File::Copy;
use File::Basename;
use Parallel::ForkManager;

use Getopt::Long;

my $filter='*';
my $in_dir;
my $out_dir;
my $help;
my $parallel=0;
GetOptions(
    "in_dir=s"=>\$in_dir,
    "out_dir=s"=>\$out_dir,
    "filter=s" => \$filter,
    "parallel=i"=>\$parallel,
    "help"           => \$help
    );

my $usage = 'run_raxml.pl [-h -f <file input filter>] -i <input_directory> -o <output_directory>' . "\n\n"; 

if($help || !defined($in_dir) || !defined($out_dir)){
    print $usage;
    exit;
}

unless(-d $in_dir){
    die "$in_dir is not a directory";

}

my $pm = new Parallel::ForkManager($parallel);

my $scratch_name= basename($in_dir);

my $name='tree';
#my $node_scratch = '/scratch/mlangill';
my $node_scratch = '~/tmp/picrust';

system("mkdir -p $node_scratch");

my $scratch_dir= $node_scratch.'/'.$scratch_name;

system("cp -r $in_dir $node_scratch");

my @input_files=glob($scratch_dir.'/'.$filter);


foreach my $input_file(@input_files){
    my $pid = $pm->start and next; 

    my ( $job_name, $input_dir) = fileparse($input_file, qr/\Q.aln\E/);
    my $output_name =  "$job_name.$name";
    my $output_file = $out_dir.'/'.$output_name;

    #Don't run job if output file already exists
    #This allows us to re-run the script to finish failed jobs only
    next if -e $output_file;

    my $local_output_file= "$scratch_dir/$output_name";

    #run the commnad
    my $command = "raxmlHPC -m GTRGAMMA -w $out_dir -s $input_file -n $output_name >/dev/null";
    system($command);

    #move the file
    #copy($local_output_file,$output_file) or die "Copy failed: $!";
    $pm->finish;
}
$pm->wait_all_children;

#remove files on node storage
system("rm -r $scratch_dir");
