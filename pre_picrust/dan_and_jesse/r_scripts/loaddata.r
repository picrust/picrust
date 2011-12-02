# LOAD DATA
tree <- read.tree('../validation/ko/data/16S_seqs_nuc_vs_GreenGenesCore_May_5.fasta_aligned_pfiltered_colons_removed.tree.noparalog')

#~     # test set with only 5 genes
#~     tipstates <- read.table('../validation/ko_small_subset/data/kegg_genome_data_April_15_2010_no_ambig_test.fna_vector.txt', row.names=1, header=F,sep='\t')

# all 4762 genes
#tipstates <- read.table('../validation/ko/data/kegg_genome_data_April_15_2010.fna_vector.txt', row.names=1, header=F,sep='\t')

# 13,000 genes
#tipstates <- read.table('../validation/full_ko_set/ko_by_genome_dump_Mar_15_2011_counts_True.txt', row.names=1, header=T)
tipstates <- read.table('../validation/full_ko_set/ko_by_genome_dump_Mar_15_2011_counts_False.txt', row.names=1, header=T)
print(rownames(tipstates)[1:10])

for(label in tree$tip.label){if(sum(rownames(tipstates)==label)==0) tree = drop.tip(tree,label)}
tipstates <- tipstates[tree$tip.label,]
tree <- multi2di(tree)
for(i in 1:length(tree$node.label)){if(tree$node.label[i]==""){tree$node.label[i] <- sprintf('inode%d',i)}}


# fix 0-length edges in tree
tree$edge.length[tree$edge.length==0] <- 1e-10
dists <- cophenetic.phylo(x=tree)
