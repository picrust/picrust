#!/usr/bin/env Rscript

Args <- commandArgs(TRUE)

dat<-read.delim(Args[1],header=FALSE);
names(dat)<-c("taxon", "values")

png("taxon_accuracy.png",width=2400,height=1200)
#png("taxon_accuracy.png")

par(mar=c(16,4,1,1))
boxplot(values ~ taxon, data=dat,las=2,ylab="Accuracy",col="bisque")
mtext(1, text=Args[2],line=10,cex.lab=18)

graphics.off()
