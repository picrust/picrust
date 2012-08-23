#!/usr/bin/env Rscript

library(vioplot)

random_precision<-read.delim("results/pfam_random_precision_accuracy.txt",header=FALSE)
random_recall<-read.delim("results/pfam_random_recall_accuracy.txt",header=FALSE)

nearest_recall<-read.delim("results/pfam_neighbour_recall_accuracy.txt",header=FALSE)
nearest_precision<-read.delim("results/pfam_neighbour_precision_accuracy.txt",header=FALSE);

ace_recall<-read.delim("results/pfam_pic_recall_accuracy.txt",header=FALSE)
ace_precision<-read.delim("results/pfam_pic_precision_accuracy.txt",header=FALSE);

wagner_recall<-read.delim("results/pfam_wagner_recall_accuracy.txt",header=FALSE)
wagner_precision<-read.delim("results/pfam_wagner_precision_accuracy.txt",header=FALSE);

cat("\tMean Median\n")
cat("PIC Recall: ",mean(ace_recall[,2]),median(ace_recall[,2]),"\n")
cat("WAGNER Recall: ",mean(wagner_recall[,2]),median(wagner_recall[,2]),"\n")
cat("NN Recall: ",mean(nearest_recall[,2]),median(nearest_recall[,2]),"\n")
cat("PIC Precision: ",mean(ace_precision[,2]),median(ace_precision[,2]),"\n")
cat("Wagner Precision: ",mean(wagner_precision[,2]),median(wagner_precision[,2]),"\n")
cat("NN Precision: ",mean(nearest_precision[,2]),median(nearest_precision[,2]),"\n")

#png("random_vs_nearest_neighbour_leave_one_out.png",width=600,height=400)
png("figures/all_methods_vs_pfam_simple_leave_one_out.png")

#vioplot(random_precision[,2],random_recall[,2],nearest_precision[,2],nearest_recall[,2],ace_precision[,2],ace_recall[,2],ylim=c(0,1),names=c("Random Precision","Random Recall","Nearest Neighbour Precision","Nearest Neighbour Recall","ACE Precision","ACE Recall"),sm_density=FALSE)

vioplot(random_precision[,2],nearest_precision[,2],ace_precision[,2],wagner_precision[,2],ylim=c(0,1),names=c("Random","Nearest Neighbour","PIC","WAGNER"),sm_density=FALSE)

title(ylab="Accuracy")

graphics.off()
