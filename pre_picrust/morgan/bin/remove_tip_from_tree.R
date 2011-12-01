#!/usr/bin/env Rscript
library(ape)
Args <- commandArgs(TRUE)

# read the tree
tree <- read.tree(Args[1])
tip_to_remove <-Args[2]
new_tree_name <-Args[3]

new_tree<-drop.tip(tree,tip_to_remove)

write.tree(new_tree,file=new_tree_name)
