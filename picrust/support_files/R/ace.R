#!/usr/bin/env Rscript
#./find_ace.R  <input: reference_tree_with_internal_nodes_labelled> <input: functional_matrix_file> <ace_method> <output: count_table> <output: prob_table>

library(ape)
Args <- commandArgs(TRUE)

#load in the tree
tree<-read.tree(Args[1])

#load in the trait table
data <-read.table(Args[2],check.names=FALSE,row.names=1,header=TRUE)

#order the trait table to match the tree tip labels
#Note: Can't just reorder with simple "data_ordered <- data[tree$tip.label,]" since this returns a vector (with no row or column names) ONLY when there is a single column in the trait table
data_ordered<-as.data.frame(data[tree$tip.label,])
rownames(data_ordered)<-tree$tip.label
names(data_ordered)<-names(data)

#do the actual ace reconsructions
reconstructions<-apply(data_ordered,2,ace,tree, type="continuous",method=Args[3])

#pull out only the ace node predictions
just_ace<-lapply(1:length(reconstructions),function(x) reconstructions[[x]]$ace)
names(just_ace)<-names(reconstructions)

#reformat the list into a matrix
just_ace_matrix<-do.call(cbind,just_ace)

#relabel the node names (ones created internally by ape) with the actual node labels in the tree
just_ace_matrix<-cbind(tree$node.label,just_ace_matrix)

#give a simple header label for the internal nodes
colnames(just_ace_matrix)[1]<-'nodes'

#Convert to a data frame
out_matrix<-data.frame(just_ace_matrix,check.names=FALSE)

#write to file
write.table(out_matrix,file=Args[4],row.names=FALSE,quote=FALSE, sep="\t")
 
