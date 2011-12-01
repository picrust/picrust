#!/usr/bin/perl
use warnings;
use strict;
use autodie;
use JSON;
use Bio::TreeIO;
use IO::String;
use Log::Log4perl;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

my $abs_dir=dirname(abs_path($0)).'/';
$abs_dir=dirname($abs_dir).'/';

my $method;
my $func;
GetOptions("method=s" =>\$method, "func=s"=>\$func);

# Set up the logger             
my $logger_cfg = $abs_dir."logger.conf";
Log::Log4perl::init($logger_cfg);
my $logger = Log::Log4perl->get_logger;

my $query_id_name=$ARGV[0];
my $ref_id_name=$ARGV[1];

my $usage ="ace.pl --method [ace|aceML] --func [pfam|subsystem|role|EC] <query_id_name>";
unless (defined($method) && ($func)){
    print $usage;
    exit;
}


my $query_id_dir =$abs_dir.'tmp/'.$query_id_name.'/';
my $neighbour_output_file=$query_id_dir.$func.'_'.$method.'_predictions.txt';
open(my $OUT,'>',$neighbour_output_file);

$logger->info("Running $0 and will store output in $neighbour_output_file");

#get pfam counts for reference genomes
my $genome_counts_file=$abs_dir.'data/'.$func.'_vs_gp_id.txt';
my $genome_counts = load_R_matrix_as_hash_of_hash($genome_counts_file);

#get ace counts for reference genomes
my $ref_dir = $abs_dir.'ref_trees/'.$ref_id_name;
my $ace_count_file=$ref_dir.'/'.$func.'_'.$method.'_counts.txt.gz';
my $ace_counts = load_R_matrix_as_hash_of_hash($ace_count_file);

#Get a list of all PFAM ids
my %pfam_id_uniq;
foreach my $gp_id(keys %{$genome_counts}){
    foreach (keys %{$genome_counts->{$gp_id}}){
        $pfam_id_uniq{$_}=1;
    }
}
my @pfam_ids=keys %pfam_id_uniq;

#load the reference tree into bioperl
my $ref_tree_file=$ref_dir.'/RAxML_result.16s_with_node_labels';
my $treeio = Bio::TreeIO->new(-file=>$ref_tree_file,-format=>'newick');
my $ref_tree=$treeio->next_tree();


#parse tree from json output into a string
my $pplacer_json_file=$query_id_dir.'pynast_trimmed_alignment.jplace';
open(my $JSON,'<',$pplacer_json_file);
my $json_text = do {local $/; <$JSON> };
my $pplacer=decode_json $json_text;
my $pplacer_tree=$pplacer->{'tree'};

#Get pplacer's highest branch label (always the last branch id([xxx]) in the tree)
my $max_branch_label;
if($pplacer_tree =~/.*\[(\d+)\]\;/){
    $max_branch_label = $1;
}else{
    $logger->fatal("can't get biggest branch label from: $pplacer_tree");
    die;
}

#parse the branch labels and get ancestor and descendant for each branch
#store information for each branch label (e.g. each read) in an array of hashes
my @nearest_node;
for my $branch_label (0..$max_branch_label){
    my $descendant_label;
    if($pplacer_tree =~ /(\w+)\:[+-]?\ *(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?\[($branch_label)\]/){
        $descendant_label=$1;
    }else{
	$logger->error("Can't get descendant node for branch: $branch_label");
	die;
    }
    my $descendant = $ref_tree->find_node(-id => $descendant_label);
   
    if(defined($descendant)){
	my $ancestor = $descendant->ancestor();
	my $branch_length = $descendant->branch_length();
	unless(defined($ancestor)){
	    $logger->warn("Can't find ancestor for node: $descendant_label and branch: $branch_label");
	    $nearest_node[$branch_label]={'descendant'=>$descendant->id(),'is_leaf'=>$descendant->is_Leaf(),'good_branch'=>1};
	}else{
	    my $ancestor_id = $ancestor->id();
	    $nearest_node[$branch_label]={'descendant'=>$descendant->id(), 'ancestor'=>$ancestor_id,'branch_length'=>$branch_length,'is_leaf'=>$descendant->is_Leaf(),'good_branch'=>2};
	}
    }else{
	$logger->warn("Can't find descendant for node: $descendant_label and branch: $branch_label");
	$nearest_node[$branch_label]={'good_branch'=>0};
    }
}


###Do the actual pfam predictions
my %pfam_predictions; #2 key hash (1st key: read id, 2nd key pfam id) value: predicted pfam count
foreach my $placement (@{$pplacer->{'placements'}}){

    my $best_place;
    foreach my $place (@{$placement->{'p'}}){
	if (!defined($best_place) || $place->[2] < $best_place->[2]){
	    next if $nearest_node[$place->[0]]{'good_branch'} ==0;
	    $best_place=$place;
	}
    }

    unless($best_place){
	my $read_str=join(',',@{$placement->{'n'}});
	$logger->error("Couldn't place reads: $read_str");
	next;
    }
    my $branch_label=$best_place->[0];
    my ($distal_length,$branch_length,$ancestor_weight,$descendant_weight);
   if($nearest_node[$branch_label]{'good_branch'}==2){
       $distal_length=$best_place->[3];
       $branch_length=$nearest_node[$branch_label]{'branch_length'};
       $ancestor_weight=$distal_length/$branch_length;
       $descendant_weight=($branch_length -$distal_length)/$branch_length;
   }else{
       $logger->warn("No ancestor information for branch: $branch_label . Using count directly from descendant node."); 
   }

    foreach my $pfam_id (@pfam_ids){
	my $descendant_count;
	if($nearest_node[$branch_label]{'is_leaf'}){
	    #use counts from genome since the descendant is a real genome
	    $descendant_count=$genome_counts->{$nearest_node[$branch_label]{'descendant'}}{$pfam_id};
	}else{
	    #use ace counts since the descendant is an internal node
	    $descendant_count=$ace_counts->{$nearest_node[$branch_label]{'descendant'}}{$pfam_id};
	}
	my $weighted_count;
	if($nearest_node[$branch_label]{'good_branch'}==2){
	    my $ancestor_count= $ace_counts->{$nearest_node[$branch_label]{'ancestor'}}{$pfam_id};
	    $weighted_count=($ancestor_weight*$ancestor_count) + ($descendant_weight*$descendant_count);
	}else{
	    $weighted_count=$descendant_count;
	}
	#round to the nearest integer
	$weighted_count=int($weighted_count+0.5);
	
	#can be more than one read mapping to the exact same branch so store the counts for all of these reads
	foreach my $read_id (@{$placement->{'n'}}){
	    $pfam_predictions{$read_id}{$pfam_id}=$weighted_count;
	}
    }
}


###output the results###
my @read_ids = sort(keys %pfam_predictions);
#matrix delimiter
my $delim="\t";

#output the header sample line
my $header = join($delim,@read_ids);
print $OUT $header,"\n";

#output each row of pfams
foreach my $pfam_id (@pfam_ids){

    my @row_counts=map {$pfam_predictions{$_}{$pfam_id}}@read_ids;
    my $row = join($delim,$pfam_id,@row_counts);
    print $OUT $row,"\n";
}


sub load_R_matrix_as_hash_of_hash{
    my $file=shift;
    my $matrix;
    my $MATRIX;
    if($file =~ /\.gz$/){
	open($MATRIX,"zcat $file|")  || $logger->fatal("Couldn't open matrix file: $file");
    }else{
	open($MATRIX,'<',$file)  || $logger->fatal("Couldn't open matrix file: $file");
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
