#!/usr/bin/env Rscript


minNeighbor<- function(tnode, pairdists,ttt){
	min_tnode <-which(pairdists[tnode,pairdists[tnode,]>0]==min(pairdists[tnode,pairdists[tnode,]>0]))
	#if there are two identical ref sequences then there can be 2 min_tnodes. Arbitralily choose the first one.
	min_tnode<-min_tnode[1]
	tnode_distance<-pairdists[tnode,min_tnode]
	read_name<-ttt$tip.label[tnode]
	ref_name<-names(min_tnode)
	return (c(read_name,ref_name,tnode_distance))
	#cat(read_name,ref_name,tnode_distance,"\n")
}


library(ape)
Args <- commandArgs(TRUE)

# read the tree
ttt <- read.tree(Args[1])
read_ids <-read.table(Args[2],check.names=FALSE)
read_ids<-as.vector(read_ids$V1)
#find the target seq
tnodes <- which(ttt$tip.label %in% read_ids)

# compute all pairwise distances to other nodes
# find the closest node to our target
pairdists<-cophenetic(ttt)

#set all pairs of distances between reads to something arbtralily large
pairdists[tnodes,tnodes]<-max(pairdists)+1
	
output_table<-sapply(tnodes, minNeighbor, pairdists,ttt)
better_format_output<-t(output_table)
write.table(better_format_output,file="",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


