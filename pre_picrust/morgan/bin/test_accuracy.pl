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

my $query_id_name=$ARGV[0];

my $usage ="test_accuracy.pl --method [random|ace|neighbour|aceML|wagner] --func [pfam|subsystem|role|EC] <query_id_name>";
unless (defined($method) && ($func)){
    print $usage;
    exit;
}


#get pfam counts for reference genomes
my $genome_counts_file=$abs_dir.'data/'.$func.'_vs_gp_id.txt';
my $genome_counts = load_R_matrix_as_hash_of_hash($genome_counts_file);


my $query_id_dir =$abs_dir.'tmp/'.$query_id_name.'/';

#choices are 'neighbour' or 'random' or 'ace'
my $prediction_type=$method;

die unless $prediction_type;

my $pfam_predictions_file = $query_id_dir.$func.'_'.$prediction_type.'_predictions.txt';
my $predicted_pfams=load_R_matrix_as_hash_of_hash($pfam_predictions_file);


my %pfam_id_uniq;
foreach my $gp_id(keys %{$genome_counts}){
    foreach (keys %{$genome_counts->{$gp_id}}){
        $pfam_id_uniq{$_}=1;
    }
}
my @pfam_ids=keys %pfam_id_uniq;

#fix bug with pfam_count.storable (some genomes don't have right gp_id)
#just remove them here (those genomes with no pfam predictions)
foreach my $gp_id (keys %{$predicted_pfams}){
    delete $predicted_pfams->{$gp_id} unless sum(values(%{$predicted_pfams->{$gp_id}}));
}

calc_proportional_accuracy();

#calculate average precision and recall for each genome (precision and recall are calculated for every pfam count)
#pfam counts where known AND predicted are both 0 are not included in calculation
#Precision = TP/(TP+FP)
#Recall = TP/(TP+FN)
sub calc_proportional_accuracy{
    my $recall_file = $query_id_dir.$func.'_'.$prediction_type.'_recall_accuracy.txt';
    my $precision_file = $query_id_dir.$func.'_'.$prediction_type.'_precision_accuracy.txt';
    
    open(my $RECALL,'>',$recall_file);
    open(my $PRECISION, '>',$precision_file);

    foreach my $gp_id(keys(%$predicted_pfams)){
	next unless (exists($genome_counts->{$gp_id}) && exists($predicted_pfams->{$gp_id}));

	my @precisions;
	my @recalls;
	foreach my $pfam_id(@pfam_ids){
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
	    push @precisions,precision(\%accuracy);
	    push @recalls,recall(\%accuracy);

	}
	my $avg_precision=sum(@precisions)/scalar(@precisions);
	my $avg_recall = sum(@recalls)/scalar(@recalls);
	print $PRECISION join("\t",$gp_id,sprintf("%.3f",$avg_precision)),"\n";
	print $RECALL join("\t",$gp_id,sprintf("%.3f",$avg_recall)),"\n";
    }
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




