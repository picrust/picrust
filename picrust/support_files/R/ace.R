#!/usr/bin/env Rscript
#./find_ace.R  <input: reference_tree_with_internal_nodes_labelled> <input: functional_matrix_file> <ace_method> <output: count_table> <output: prob_table>

library(ape)
Args <- commandArgs(TRUE)

#load in the tree
tree<-read.tree(Args[1])

#load in the trait table
data <-read.table(Args[2],check.names=FALSE,row.names=1,header=TRUE)

#sort the data matrix so it matches the tree (this is needed before running ace())
data <- data[tree$tip.label, ]

#do the actual ace reconsructions
reconstructions<-apply(data,2,ace,tree, type="continuous",method=Args[3])

#pull out only the ace node predictions
just_ace<-lapply(1:length(reconstructions),function(x) reconstructions[[x]]$ace)
names(just_ace)<-names(reconstructions)

#reformat the list into a matrix
just_ace_matrix<-do.call(rbind,just_ace)

#relabel the node names (ones created internally by ape) with the actual node labels in the tree
colnames(just_ace_matrix)<-tree$node.label

just_ace_matrix<-t(just_ace_matrix)

nodes<-row.names(just_ace_matrix)

out_matrix<-data.frame(nodes,just_ace_matrix,check.names=FALSE)
#write the matrix to file
write.table(out_matrix,file=Args[4],row.names=FALSE,quote=FALSE, sep="\t")
 