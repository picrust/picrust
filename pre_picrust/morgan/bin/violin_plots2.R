#!/usr/bin/env Rscript

library(vioplot)

#load in files
pfam_random_precision<-read.delim("results/pfam_random_precision_accuracy.txt",header=FALSE)
pfam_random_recall<-read.delim("results/pfam_random_recall_accuracy.txt",header=FALSE)
pfam_nearest_recall<-read.delim("results/pfam_neighbour_recall_accuracy.txt",header=FALSE)
pfam_nearest_precision<-read.delim("results/pfam_neighbour_precision_accuracy.txt",header=FALSE);
pfam_pic_recall<-read.delim("results/pfam_pic_recall_accuracy.txt",header=FALSE)
pfam_pic_precision<-read.delim("results/pfam_pic_precision_accuracy.txt",header=FALSE);

role_random_precision<-read.delim("results/role_random_precision_accuracy.txt",header=FALSE)
role_random_recall<-read.delim("results/role_random_recall_accuracy.txt",header=FALSE)
role_nearest_recall<-read.delim("results/role_neighbour_recall_accuracy.txt",header=FALSE)
role_nearest_precision<-read.delim("results/role_neighbour_precision_accuracy.txt",header=FALSE);
role_pic_recall<-read.delim("results/role_pic_recall_accuracy.txt",header=FALSE)
role_pic_precision<-read.delim("results/role_pic_precision_accuracy.txt",header=FALSE);

subsystem_random_precision<-read.delim("results/subsystem_random_precision_accuracy.txt",header=FALSE)
subsystem_random_recall<-read.delim("results/subsystem_random_recall_accuracy.txt",header=FALSE)
subsystem_nearest_recall<-read.delim("results/subsystem_neighbour_recall_accuracy.txt",header=FALSE)
subsystem_nearest_precision<-read.delim("results/subsystem_neighbour_precision_accuracy.txt",header=FALSE);
subsystem_pic_recall<-read.delim("results/subsystem_pic_recall_accuracy.txt",header=FALSE)
subsystem_pic_precision<-read.delim("results/subsystem_pic_precision_accuracy.txt",header=FALSE);

EC_random_precision<-read.delim("results/EC_random_precision_accuracy.txt",header=FALSE)
EC_random_recall<-read.delim("results/EC_random_recall_accuracy.txt",header=FALSE)
EC_nearest_recall<-read.delim("results/EC_neighbour_recall_accuracy.txt",header=FALSE)
EC_nearest_precision<-read.delim("results/EC_neighbour_precision_accuracy.txt",header=FALSE);
EC_pic_recall<-read.delim("results/EC_pic_recall_accuracy.txt",header=FALSE)
EC_pic_precision<-read.delim("results/EC_pic_precision_accuracy.txt",header=FALSE);

#cat("\tMean Median\n")
#cat("ACE Recall: ",mean(ace_recall[,2]),median(ace_recall[,2]),"\n")
#cat("Near Recall: ",mean(nearest_recall[,2]),median(nearest_recall[,2]),"\n")
#cat("ACE Precision: ",mean(ace_precision[,2]),median(ace_precision[,2]),"\n")
#cat("Near Precision: ",mean(nearest_precision[,2]),median(nearest_precision[,2]),"\n")

png("figures/leave_one_out_accuracy_of_all_methods_on_all_funcs.png",width=2000,height=1200)
plot.new()
plot.window(xlim=c(0.2,24.2),ylim=c(0,1))
par(mar = c(7, 4, 4, 2) + 0.1)
vioplot(pfam_random_precision[,2],at=1,add=TRUE,sm_density=FALSE)
vioplot(pfam_random_recall[,2],at=2,add=TRUE,sm_density=FALSE)
vioplot(pfam_nearest_precision[,2],at=3,add=TRUE,sm_density=FALSE)
vioplot(pfam_nearest_recall[,2],at=4,add=TRUE,sm_density=FALSE)
vioplot(pfam_pic_precision[,2],at=5,add=TRUE,sm_density=FALSE)
vioplot(pfam_pic_recall[,2],at=6,add=TRUE,sm_density=FALSE)

vioplot(role_random_precision[,2],at=7,add=TRUE,sm_density=FALSE)
vioplot(role_random_recall[,2],at=8,add=TRUE,sm_density=FALSE)
vioplot(role_nearest_precision[,2],at=9,add=TRUE,sm_density=FALSE)
vioplot(role_nearest_recall[,2],at=10,add=TRUE,sm_density=FALSE)
vioplot(role_pic_precision[,2],at=11,add=TRUE,sm_density=FALSE)
vioplot(role_pic_recall[,2],at=12,add=TRUE,sm_density=FALSE)

vioplot(subsystem_random_precision[,2],at=13,add=TRUE,sm_density=FALSE)
vioplot(subsystem_random_recall[,2],at=14,add=TRUE,sm_density=FALSE)
vioplot(subsystem_nearest_precision[,2],at=15,add=TRUE,sm_density=FALSE)
vioplot(subsystem_nearest_recall[,2],at=16,add=TRUE,sm_density=FALSE)
vioplot(subsystem_pic_precision[,2],at=17,add=TRUE,sm_density=FALSE)
vioplot(subsystem_pic_recall[,2],at=18,add=TRUE,sm_density=FALSE)

vioplot(EC_random_precision[,2],at=19,add=TRUE,sm_density=FALSE)
vioplot(EC_random_recall[,2],at=20,add=TRUE,sm_density=FALSE)
vioplot(EC_nearest_precision[,2],at=21,add=TRUE,sm_density=FALSE)
vioplot(EC_nearest_recall[,2],at=22,add=TRUE,sm_density=FALSE)
vioplot(EC_pic_precision[,2],at=23,add=TRUE,sm_density=FALSE)
vioplot(EC_pic_recall[,2],at=24,add=TRUE,sm_density=FALSE)

#Add ticks and labels on y axis
axis(2)


x_labels<-c(
"PFAM R P",
"PFAM R R",
"PFAM NN P",
"PFAM NN R",
"PFAM PIC P",
"PFAM PIC R",
"Role R P",
"Role R R",
"Role NN P",
"Role NN R",
"Role PIC P",
"Role PIC R",
"SS R P",
"SS R R",
"SS NN P",
"SS NN R",
"SS PIC P",
"SS PIC R",
"EC R P",
"EC R R",
"EC NN P",
"EC NN R",
"EC PIC P",
"EC PIC R"
)


#add ticks only to x axis
axis(1,1:24,x_labels,las=2)
#axis(1,labels=FALSE)
 
## Plot x axis labels at default tick marks
#text(1:6, par("usr")[3] - 0.25, srt = 45, adj = 1, labels = x_labels, xpd = TRUE)


graphics.off()
