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

my $query_id_name=$ARGV[0];

my $query_id_dir =$abs_dir.'tmp/'.$query_id_name.'/';
my $neighbour_output_file=$query_id_dir.$func.'_neighbour_predictions.txt';
open(my $OUT,'>',$neighbour_output_file);

#Find the nearest neighbour for each read using find_neighbour.R
my $read_ids_file=$query_id_dir.'read_ids.txt';
my $pplacer_tree_file = $query_id_dir.'pplacer.tog.tre';
my $find_neighbour_file = $query_id_dir.'find_neighbour.txt';
my $find_neighbour_cmd=$abs_dir."bin/find_neighbour.R $pplacer_tree_file $read_ids_file > $find_neighbour_file";
system($find_neighbour_cmd);

#get pfam counts for reference genomes
my $genome_counts_file=$abs_dir.'data/'.$func.'_vs_gp_id.txt';
my $genome_counts = load_R_matrix_as_hash_of_hash($genome_counts_file);


open (my $IN, $find_neighbour_file);
my %ref_map;
while (<$IN>){
    chomp;
    my ($read_id,$ref,$dist)=split;

    $ref_map{$read_id}=$ref;
}

#Get a list of all PFAM ids
my %pfam_id_uniq;
foreach my $gp_id(keys %{$genome_counts}){
    foreach (keys %{$genome_counts->{$gp_id}}){
        $pfam_id_uniq{$_}=1;
    }
}
my @pfam_ids=keys %pfam_id_uniq;


my @read_ids = sort(keys %ref_map);


#matrix delimiter
my $delim="\t";

#output the header sample line
my $header = join($delim,@read_ids);
print $OUT $header,"\n";

#output each row of pfams
foreach my $pfam_id (@pfam_ids){
    my @row_counts;
    foreach my $read_id(@read_ids){
	my $gp_id=$ref_map{$read_id};
        if(exists($genome_counts->{$gp_id}{$pfam_id})){
	    push @row_counts,$genome_counts->{$gp_id}{$pfam_id};
	}else{
	    push @row_counts,0;
	}

    }
    my $row = join($delim,$pfam_id,@row_counts);
    print $OUT $row,"\n";
}


sub load_R_matrix_as_hash_of_hash{
    my $file=shift;
    my $matrix;
    open(my $MATRIX,'<',$file);
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
