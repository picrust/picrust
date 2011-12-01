#!/usr/bin/env Rscript
#./find_ace.R <output: file.gz> <input: functional_matrix_file> <input: reference_tree> <ace_method>

library(ape)
library(geiger)
Args <- commandArgs(TRUE)

out_file_FH<-gzfile(Args[1],"w")

tree<-read.tree(Args[3])

#give the tree node labels (same as those used by Phylo internally)
tree$node.label<-Ntip(tree)+(1:(Ntip(tree)-1))

new_tree_file<-paste(Args[3],'_with_node_labels',sep="")
#and write it back to file
write.tree(tree,file=new_tree_file)

tree<-multi2di(tree,random=FALSE)


#load in the pfam data
data <-read.table(Args[2],check.names=FALSE)
data<-t(data)

#check for name overlaps in tree and data structures
overlap<-name.check(tree,data)

#get rid of taxa that we don't have data for
#(assume we don't have to do this)
#tree<-drop.tip(tree,overlap$Tree.not.data)

#get rid of data that are not in our reference tree
data<-data[!(rownames(data) %in% overlap$Data.not.tree),]

#sort the data matrix so it matches the tree (this is needed before running ace())
data <- data[tree$tip.label, ]

#do the actual ace reconsructions
pic_reconstructions<-apply(data,2,ace,tree, type="continuous",method=Args[4],CI=FALSE)

#pull out only the ace node predictions
just_ace<-lapply(1:length(pic_reconstructions),function(x) pic_reconstructions[[x]]$ace)
names(just_ace)<-names(pic_reconstructions)

#reformat the list into a matrix
just_ace_matrix<-do.call(rbind,just_ace)

#write the matrix to file
write.table(just_ace_matrix,file=out_file_FH,quote=FALSE, sep="\t")
 
close(out_file_FH)