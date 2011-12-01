#!/usr/bin/perl
use warnings;
use strict;
use autodie;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

my $abs_dir=dirname(abs_path($0)).'/';
$abs_dir=dirname($abs_dir).'/';


my $func;
GetOptions("func=s"=>\$func);

my $output_file=$abs_dir.'data/'.$func.'_metadata.txt';
open(my $OUT,'>',$output_file);

#get pfam counts for reference genomes
my $genome_counts_file=$abs_dir.'data/'.$func.'_vs_gp_id.txt';
my $genome_counts = load_R_matrix_as_hash_of_hash($genome_counts_file);

#Get a list of all PFAM ids
my %pfam;
foreach my $gp_id(keys %{$genome_counts}){
    foreach my $pfam_id(keys %{$genome_counts->{$gp_id}}){
	if($genome_counts->{$gp_id}{$pfam_id}>0){
	    if(exists($pfam{$pfam_id}{'abundance'})){
		$pfam{$pfam_id}{'abundance'}+=$genome_counts->{$gp_id}{$pfam_id};
	    }else{
		$pfam{$pfam_id}{'abundance'}=$genome_counts->{$gp_id}{$pfam_id};
	    }
	    $pfam{$pfam_id}{coverage}++;
	    if(!exists($pfam{$pfam_id}{max}) || $pfam{$pfam_id}{max}<$genome_counts->{$gp_id}{$pfam_id}){
		$pfam{$pfam_id}{max}=$genome_counts->{$gp_id}{$pfam_id};
	    }
	}
    }
}

my $delim="\t";
my $header="description\tabundance\tcoverage\tmax\n";
print $OUT $header;

my $descriptions_file=$abs_dir.'data/'.$func.'_descriptions.txt';
open(my $NAMES,'<',$descriptions_file);

while(<$NAMES>){
    chomp;
    my($id,$desc)=split(/\t/,$_);
    next unless exists($pfam{$id});
    my $row=join($delim,$id,$desc,$pfam{$id}{'abundance'},$pfam{$id}{coverage},$pfam{$id}{max});
    print $OUT $row,"\n";

}

sub load_R_matrix_as_hash_of_hash{
    my $file=shift;
    my $matrix;
    my $MATRIX;
    if($file =~ /\.gz$/){
	open($MATRIX,"zcat $file|")  || die("Couldn't open matrix file: $file");
    }else{
	open($MATRIX,'<',$file)  || die("Couldn't open matrix file: $file");
    }
    my $header = <$MATRIX>;
    chomp($header);
    my (@col_names)=split(/\s+/,$header);
    while(<$MATRIX>){
	chomp;
	my ($row_name,@cells)=split;
	for my $index(0..$#col_names){
	    $matrix->{$col_names[$index]}{$row_name}=$cells[$index];
	}
    }
    return $matrix;
}
