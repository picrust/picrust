#!/usr/bin/env perl
#perldoc make_predictions.pl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Bio::TreeIO;
use JSON;
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
GetOptions (\%opt,'func=s','method=s','ref_tree=s','query=s','Force','help') or pod2usage(2);
pod2usage(-verbose=>2) if exists $opt{'help'};

#Mandatory options
foreach my $option ('ref_tree','method','func','query'){
    pod2usage($0.': You must specify option for --'.$option) unless exists $opt{$option};
}


my $query_id_dir =$abs_dir.'tmp/'.$opt{'query'}.'/';
my $func=$opt{'func'};
my $method=$opt{'method'};
my $ref_name=$opt{'ref_tree'};

my $output_file=$query_id_dir.$func.'_'.$method.'_predictions.txt';


if(-e $output_file && ! exists $opt{'Force'}){
    $logger->info("Not going to make predictions because output file already exists. Use --Force to overwrite existing files.");
    exit;
}

$logger->info("Running $0 and will store output in $output_file");

#get functional counts for reference genomes
my $genome_counts_file=$abs_dir.'data/'.$func.'_vs_gp_id.txt';
my $genome_counts = load_R_matrix_as_hash_of_hash($genome_counts_file);

#Get a list of all PFAM ids
my %pfam_id_uniq;
foreach my $gp_id(keys %{$genome_counts}){
    foreach (keys %{$genome_counts->{$gp_id}}){
        $pfam_id_uniq{$_}=1;
    }
}
my @pfam_ids=keys %pfam_id_uniq;

if($method eq 'ML' || $method eq 'REML' || $method eq 'pic'){
    my $ace_predictions = make_ace_predictions();
    output_predictions($ace_predictions);
}elsif($method eq 'random'){
    my $random_predictions = make_random_predictions();
    output_predictions($random_predictions);
}elsif($method eq 'neighbour'){
    my $nn_predictions = make_nn_predictions();
    output_predictions($nn_predictions);
}else{
    $logger->fatal("Not valid method option: $method");
    exit(1);
}

sub make_nn_predictions{

    #Find the nearest neighbour for each read using find_neighbour.R
    my $read_ids_file=$query_id_dir.'read_ids.txt';
    my $pplacer_tree_file = $query_id_dir.'pplacer.tog.tre';
    my $find_neighbour_file = $query_id_dir.'find_neighbour.txt';
    my $find_neighbour_cmd=$abs_dir."bin/find_neighbour.R $pplacer_tree_file $read_ids_file > $find_neighbour_file";
    system($find_neighbour_cmd);


    open (my $IN, '<',$find_neighbour_file);
    my %ref_map;
    while (<$IN>) {
	chomp;
	my ($read_id,$ref,$dist)=split;

	$ref_map{$read_id}=$ref;
    }
    close($IN);

    my @read_ids = sort(keys %ref_map);

    my %nn_predictions;
    #output each row of pfams
    foreach my $pfam_id (@pfam_ids) {
	my @row_counts;
	foreach my $read_id (@read_ids) {
	    my $gp_id=$ref_map{$read_id};
	    if (exists($genome_counts->{$gp_id}{$pfam_id})) {
		$nn_predictions{$read_id}{$pfam_id}=$genome_counts->{$gp_id}{$pfam_id};
	    } else {
		$nn_predictions{$read_id}{$pfam_id}=0;
	    }

	}
    }
    return \%nn_predictions;
}

sub make_random_predictions{

    my $read_ids_file=$query_id_dir.'read_ids.txt';
    open(my $READS,'<',$read_ids_file);
    my @read_ids=<$READS>;
    chomp(@read_ids);
    close($READS);

    #For leave one out testing we want to remove(s) the genome being tested from the reference database
    foreach my $id (@read_ids) {
	if (exists $genome_counts->{$id}) {
	    delete $genome_counts->{$id};
	}
    }

    my @ref_ids=keys %$genome_counts;

    my %random_predictions;
    #output each row of pfams
    foreach my $pfam_id (@pfam_ids) {
	my @row_counts;
	foreach my $read_id (@read_ids) {
	    #get a random gp_id
	    my $gp_id = $ref_ids[ rand @ref_ids ];
	    if (exists $genome_counts->{$gp_id}{$pfam_id}) {
		$random_predictions{$read_id}{$pfam_id}=$genome_counts->{$gp_id}{$pfam_id};
	    } else {
		$random_predictions{$read_id}{$pfam_id}=0;
	    }
	}
    }
    return \%random_predictions;
}

sub make_ace_predictions{

    #get ace counts for reference genomes
    my $ref_dir = $abs_dir.'ref_trees/'.$ref_name;
    my $ace_count_file=$ref_dir.'/'.$func.'_'.$method.'_counts.txt.gz';
    my $ace_counts = load_R_matrix_as_hash_of_hash($ace_count_file);


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
    if ($pplacer_tree =~/.*\[(\d+)\]\;/) {
	$max_branch_label = $1;
    } else {
	$logger->fatal("can't get biggest branch label from: $pplacer_tree");
	die;
    }

    #parse the branch labels and get ancestor and descendant for each branch
    #store information for each branch label (e.g. each read) in an array of hashes
    my @nearest_node;
    for my $branch_label (0..$max_branch_label) {
	my $descendant_label;
	if ($pplacer_tree =~ /(\w+)\:[+-]?\ *(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?\[($branch_label)\]/) {
	    $descendant_label=$1;
	} else {
	    $logger->error("Can't get descendant node for branch: $branch_label");
	    die;
	}
	my $descendant = $ref_tree->find_node(-id => $descendant_label);
   
	if (defined($descendant)) {
	    my $ancestor = $descendant->ancestor();
	    my $branch_length = $descendant->branch_length();
	    unless(defined($ancestor)){
		$logger->warn("Can't find ancestor for node: $descendant_label and branch: $branch_label");
		$nearest_node[$branch_label]={'descendant'=>$descendant->id(),'is_leaf'=>$descendant->is_Leaf(),'good_branch'=>1};
	    } else {
		my $ancestor_id = $ancestor->id();
		$nearest_node[$branch_label]={'descendant'=>$descendant->id(), 'ancestor'=>$ancestor_id,'branch_length'=>$branch_length,'is_leaf'=>$descendant->is_Leaf(),'good_branch'=>2};
	    }
	} else {
	    $logger->warn("Can't find descendant for node: $descendant_label and branch: $branch_label");
	    $nearest_node[$branch_label]={'good_branch'=>0};
	}
    }


    ###Do the actual pfam predictions
    my %pfam_predictions; #2 key hash (1st key: read id, 2nd key pfam id) value: predicted pfam count
    foreach my $placement (@{$pplacer->{'placements'}}) {

	my $best_place;
	foreach my $place (@{$placement->{'p'}}) {
	    if (!defined($best_place) || $place->[2] < $best_place->[2]) {
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
	if ($nearest_node[$branch_label]{'good_branch'}==2) {
	    $distal_length=$best_place->[3];
	    $branch_length=$nearest_node[$branch_label]{'branch_length'};
	    $ancestor_weight=$distal_length/$branch_length;
	    $descendant_weight=($branch_length -$distal_length)/$branch_length;
	} else {
	    $logger->warn("No ancestor information for branch: $branch_label . Using count directly from descendant node."); 
	}

	foreach my $pfam_id (@pfam_ids) {
	    my $descendant_count;
	    if ($nearest_node[$branch_label]{'is_leaf'}) {
		#use counts from genome since the descendant is a real genome
		$descendant_count=$genome_counts->{$nearest_node[$branch_label]{'descendant'}}{$pfam_id};
	    } else {
		#use ace counts since the descendant is an internal node
		$descendant_count=$ace_counts->{$nearest_node[$branch_label]{'descendant'}}{$pfam_id};
	    }
	    my $weighted_count;
	    if ($nearest_node[$branch_label]{'good_branch'}==2) {
		my $ancestor_count= $ace_counts->{$nearest_node[$branch_label]{'ancestor'}}{$pfam_id};
		$weighted_count=($ancestor_weight*$ancestor_count) + ($descendant_weight*$descendant_count);
	    } else {
		$weighted_count=$descendant_count;
	    }
	    #round to the nearest integer
	    $weighted_count=int($weighted_count+0.5);
	
	    #can be more than one read mapping to the exact same branch so store the counts for all of these reads
	    foreach my $read_id (@{$placement->{'n'}}) {
		$pfam_predictions{$read_id}{$pfam_id}=$weighted_count;
	    }
	}
    }
    return \%pfam_predictions;
}


sub output_predictions {
    my $predictions = shift;
    
    ###output the results###
    my @read_ids = sort(keys %$predictions);
    #matrix delimiter
    my $delim="\t";
    
    open(my $OUT,'>',$output_file);
    
    #output the header sample line
    my $header = join($delim,@read_ids);
    print $OUT $header,"\n";
    
    #output each row of pfams
    foreach my $pfam_id (@pfam_ids){
	
	my @row_counts=map {$predictions->{$_}{$pfam_id}}@read_ids;
	my $row = join($delim,$pfam_id,@row_counts);
	print $OUT $row,"\n";
    }
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

__END__

=head1 Name

make_predictions.pl - Makes functional metagenomic predictions.

=head1 USAGE

make_predictions.pl -f=[pfam|EC|subsystem|role] -m=[pic|ml|reml|random|neighbour] -r=<reference_tree_name> -q=<query_name> [OPTIONS]

E.g.:

make_predictions.pl -f=pfam -m=pic -r=16s_seqs_in_seed.fa -q=query.fa

=head1 OPTIONS

=over 4

=item B<-m, --method (REQUIRED)>

Defines which prediction method to run. Choices include: pic,ml,reml,random,neighbour

=item B<-f, --func (REQUIRED)>

Defines which functional traits are being used for reconstruction. Choices include: subsystem, pfam, EC, roles

=item B<-r, --ref_tree (REQUIRED)>

Specify reference tree. Note this is just the directory name within (ref_trees).

=item B<-q, --query (REQUIRED)>

Specify name of query. This is the same name as your query file used in place reads and is the directory in(tmp).

=item B<-F,--Force>

Forces all steps even if output files already exist (off by default) 

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<make_predictions.pl>  Makes functional metagenomic predictions given a reference 16S tree with reads mapped to it using a variety of methods.


=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

=head1 DATE

08-Dec-2011

=cut

