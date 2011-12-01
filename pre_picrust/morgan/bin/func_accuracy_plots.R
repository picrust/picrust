#!/usr/bin/env Rscript
Args <- commandArgs(TRUE)

func<-Args[1]

#load in data about each of the functional groups (description,coverage, abundance, max, etc.)
meta_data<-read.delim(paste("data/",func,"_metadata.txt",sep=""))

#load in the func accuracys using different methods
random_prec<-read.table(paste("results/",func,"_random_precision_func_accuracy.txt",sep=""))
near_prec<-read.table(paste("results/",func,"_neighbour_precision_func_accuracy.txt",sep=""))
ace_prec<-read.table(paste("results/",func,"_ace_precision_func_accuracy.txt",sep=""))

#merge the datasets together into a single dataframe
data<-merge(meta_data,random_prec,by=0)
data<-merge(data,near_prec,by.x="Row.names",by.y='row.names')
data<-merge(data,ace_prec,by.x="Row.names",by.y='row.names')

#make the variables (columns) of data locally accessible (e.g. 'coverage' instead of 'data$coverage')
attach(data)

#plot coverage vs accuracy
png(paste("figures/",func,"_accuracy_NN_vs_coverage.png",sep=""))
plot(coverage,neighbour_precision,xlab="Coverage",ylab="Accuracy (NN Precision)",main=paste(func," Accuracy vs Coverage",sep=""))
#calculate correlation and put it on figure
mtext(sprintf("R2 = %.3f",cor(coverage,neighbour_precision)))
graphics.off()


#if running manually you can call identify so you can click on points and get names
#identify(coverage,neighbour_precision,labels=Row.names)

#Same figure as above but with transposase PFAMs colored red
png(paste("figures/",func,"_accuracy_NN_vs_coverage_transposases.png",sep=""))
plot(coverage,neighbour_precision,xlab="Coverage",ylab="Accuracy (NN Precision)",main="Transposases show lower accuracy",pch="*",col="Grey")

####pull out "Transposase" PFAMs
transposases<-unique(c(grep("transposase",description,ignore.case=TRUE),grep("transposase",Row.names,ignore.case=TRUE)))
integrases<-grep("integrase",description,ignore.case=TRUE)
mobile_pfams<-unique(c(transposases,integrases))

phage<-unique(c(grep("phage",description,ignore.case=TRUE),grep("phage",Row.names,ignore.case=TRUE)))

points(data[phage,]$coverage,data[phage,]$neighbour_precision,col="Blue",pch="*",cex=1.5)

trans_cov<-data[mobile_pfams,]$coverage
trans_nn<-data[mobile_pfams,]$neighbour_precision

####calculate correlation and put it on figure
mtext(sprintf("R2 = %.3f",cor(trans_cov,trans_nn)))
points(trans_cov,trans_nn,col="Red",pch="*",cex=1.5)
graphics.off()

####interative plot to select particular interesting points (low accuracy, HGT pfams)
#library(iplots)
#iplot(coverage,near_precision)
##(Manually) Select points
#interesting_pfams<-iset.selected()
#names(interesting_data)[1]<-"pfam_id"
#write.table(interesting_data,"low_accuracy_pfams.txt",sep="\t",row.names=FALSE,quote=FALSE)

####plot coverage vs accuracy for ace
png(paste("figures/",func,"_accuracy_ace_vs_coverage.png",sep=""))
plot(coverage,ace_precision,xlab="Coverage",ylab="Accuracy (ACE Precision)",main=paste(func,"Accuracy vs Coverage"))
#calculate correlation and put it on figure
mtext(sprintf("R2 = %.3f",cor(coverage,ace_precision)))
graphics.off()

####plot coverage vs accuracy for random
png(paste("figures/",func,"_random_accuracy_vs_coverage.png",sep=""))
plot(coverage,random_precision,xlab="Coverage",ylab="Accuracy (Random Precision)",main=paste(func,"Accuracy vs Coverage"))
#calculate correlation and put it on figure
mtext(sprintf("R2 = %.3f",cor(coverage,random_precision)))
graphics.off()

#plot functional coverage distribution
png(paste("figures/",func,"_coverage.png",sep=""))
plot(density(coverage),xlab="Coverage (# genomes >1 functional category)",main=paste("Distribution of",func,"Coverage"))
graphics.off()

#plot functional accuracy distribution
png(paste("figures/",func,"_accuracy_distribution.png",sep=""))
plot(density(neighbour_precision),col="blue",main="Accuracy Distribution",xlab="Accuracy Across Genomes")
lines(density(ace_precision),col="red")
legend("bottom",c("NN","ACE"),col=c("blue","red"),lty=1)
graphics.off()