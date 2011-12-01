#!/usr/bin/env perl

use warnings;
use strict;
use autodie;
use List::Util qw(sum);
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

my $abs_dir=dirname(abs_path($0)).'/';
$abs_dir=dirname($abs_dir).'/';

my $method;
my $func;
GetOptions("method=s" =>\$method, "func=s"=>\$func);

#get functional counts for reference genomes
my $genome_counts_file=$abs_dir.'data/'.$func.'_vs_gp_id.txt';
my $genome_counts = load_R_matrix_as_hash_of_hash($genome_counts_file);


#choices are 'neighbour' or 'random' or 'ace'
my $prediction_type=$method;

die unless $prediction_type;

my $tmp=$abs_dir.'tmp/*/'.$func.'_'.$prediction_type.'_predictions.txt';

my @query_files=`ls -1 $tmp`;
chomp(@query_files);

my $predicted_pfams;
foreach(@query_files){
    $predicted_pfams=load_R_matrix_as_hash_of_hash($_,$predicted_pfams);
}

my %pfam_id_uniq;
foreach my $gp_id(keys %{$genome_counts}){
    foreach (keys %{$genome_counts->{$gp_id}}){
        $pfam_id_uniq{$_}=1;
    }
}
my @pfam_ids=keys %pfam_id_uniq;

#just remove them here (those genomes with no pfam predictions)
foreach my $gp_id (keys %{$predicted_pfams}){
    delete $predicted_pfams->{$gp_id} unless sum(values(%{$predicted_pfams->{$gp_id}}));
}

my %gp_id_nn;
#load in the nearest neighbour distance for each gp_id
foreach my $gp_id (keys %$predicted_pfams){
    my $nn_file= $abs_dir.'tmp/'.$gp_id.'_query.fa/find_neighbour.txt';
    open (my $NN,'<',$nn_file);
    $_=<$NN>;
    chomp;
    my($query,$nn,$dist)=split;
    $gp_id_nn{$query}=$dist;
}
calc_proportional_func_accuracy();

#calculate average precision and recall for each pfam (precision and recall are calculated for every genome...then averaged across these)
#pfam counts where known AND predicted are both 0 are not included in calculation
#Precision = TP/(TP+FP)
#Recall = TP/(TP+FN)
sub calc_proportional_func_accuracy{
    my $recall_file = $abs_dir.'results/'.$func.'_'.$prediction_type.'_recall_func_accuracy.txt';
    my $precision_file = $abs_dir.'results/'.$func.'_'.$prediction_type.'_precision_func_accuracy.txt';
    
    open(my $RECALL,'>',$recall_file);
    open(my $PRECISION, '>',$precision_file);

    my $all_recall_file = $abs_dir.'results/'.$func.'_'.$prediction_type.'_recall_func_all_accuracy.txt';
    my $all_precision_file = $abs_dir.'results/'.$func.'_'.$prediction_type.'_precision_func_all_accuracy.txt';
    
    open(my $ALL_RECALL,'>',$all_recall_file);
    open(my $ALL_PRECISION, '>',$all_precision_file);

    #print R headers
    print $PRECISION $prediction_type.'_precision',"\n";
    print $RECALL $prediction_type.'_recall',"\n";

    foreach my $pfam_id(@pfam_ids){

	my @precisions;
	my @recalls;
	
	foreach my $gp_id(keys(%$predicted_pfams)){
	    next unless (exists($genome_counts->{$gp_id}) && exists($predicted_pfams->{$gp_id}));

	    my %accuracy=(tp=>0,fp=>0,fn=>0);
	    my $known_count=0;
	    if(exists($genome_counts->{$gp_id}{$pfam_id})){
		$known_count=$genome_counts->{$gp_id}{$pfam_id};
	    }
	    my $predicted_count=0;
	    if(exists($predicted_pfams->{$gp_id}{$pfam_id})){
		$predicted_count=$predicted_pfams->{$gp_id}{$pfam_id};
	    }
	    next if ($known_count==0 && $predicted_count==0);
	    
	    if($predicted_count > $known_count){
		$accuracy{tp}+=$known_count;
		$accuracy{fp}+=$predicted_count-$known_count;
		$accuracy{fn}+=0;
	    }else{
		$accuracy{tp}+=$predicted_count;
		$accuracy{fp}+=0;
		$accuracy{fn}+=$known_count-$predicted_count;
	    }

	    #calculate precision and recall
	    my $precision =precision(\%accuracy);
	    my $recall = recall(\%accuracy);
	    push @precisions,$precision;
	    push @recalls,$recall;

	    #output these precisions and recalls
	    print $ALL_PRECISION join("\t",$gp_id,$pfam_id,sprintf("%.3f",$precision),$gp_id_nn{$gp_id}),"\n";
	    print $ALL_RECALL join("\t",$gp_id,$pfam_id,sprintf("%.3f",$recall),$gp_id_nn{$gp_id}),"\n";

	}
	my $avg_precision;
	my $avg_recall;
	if(scalar(@precisions)){
	    $avg_precision=sum(@precisions)/scalar(@precisions);
	    $avg_recall = sum(@recalls)/scalar(@recalls);
	}else{
	    #this happens when a pfam only occurs in 1 genome and it is the same genome that is being used as the query
	    #just set them manually to 0
	    $avg_precision='0';
	    $avg_recall='0';
	}
	print $PRECISION join("\t",$pfam_id,sprintf("%.3f",$avg_precision)),"\n";
	print $RECALL join("\t",$pfam_id,sprintf("%.3f",$avg_recall)),"\n";
    }
}


sub load_R_matrix_as_hash_of_hash{
    my $file=shift;
    my $matrix=shift;
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



sub precision {
    my $a=shift;
    if(($a->{tp}+$a->{fp})==0){
	return 0;
    }
    return $a->{tp}/($a->{tp}+$a->{fp});
}

sub recall {
    my $a=shift;
    if(($a->{tp}+$a->{fn})==0){
	return 0;
    }
    return $a->{tp}/($a->{tp}+$a->{fn});
}




