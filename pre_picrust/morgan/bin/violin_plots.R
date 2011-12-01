#!/usr/bin/env Rscript

library(vioplot)

random_precision<-read.delim("results/random_precision_accuracy.txt",header=FALSE)
random_recall<-read.delim("results/random_recall_accuracy.txt",header=FALSE)

nearest_recall<-read.delim("results/neighbour_recall_accuracy.txt",header=FALSE)
nearest_precision<-read.delim("results/neighbour_precision_accuracy.txt",header=FALSE);

ace_recall<-read.delim("results/ace_recall_accuracy.txt",header=FALSE)
ace_precision<-read.delim("results/ace_precision_accuracy.txt",header=FALSE);

cat("\tMean Median\n")
cat("ACE Recall: ",mean(ace_recall[,2]),median(ace_recall[,2]),"\n")
cat("Near Recall: ",mean(nearest_recall[,2]),median(nearest_recall[,2]),"\n")
cat("ACE Precision: ",mean(ace_precision[,2]),median(ace_precision[,2]),"\n")
cat("Near Precision: ",mean(nearest_precision[,2]),median(nearest_precision[,2]),"\n")

#png("random_vs_nearest_neighbour_leave_one_out.png",width=600,height=400)
png("figures/random_vs_nearest_neighbour_leave_one_out.png")

#vioplot(random_precision[,2],random_recall[,2],nearest_precision[,2],nearest_recall[,2],ace_precision[,2],ace_recall[,2],ylim=c(0,1),names=c("Random Precision","Random Recall","Nearest Neighbour Precision","Nearest Neighbour Recall","ACE Precision","ACE Recall"),sm_density=FALSE)

vioplot(random_precision[,2],nearest_precision[,2],ace_precision[,2],ylim=c(0,1),names=c("Random","Nearest Neighbour","ACE"),sm_density=FALSE)

title(ylab="Accuracy")

graphics.off()
