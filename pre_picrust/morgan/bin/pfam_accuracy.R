#!/usr/bin/env Rscript

#load in data about each of the pfams (coverage, abundance, max, etc.)
pfam<-read.table("../data/pfam_metadata.txt");

#load in more data about each of the pfams (pfam_accession, pfam_description, pfam_length)
pfam_meta<-read.delim("../data/more_pfam_metadata.txt");

 
#load in the pfam accuracys using different methods
random_prec<-read.table("results/random_precision_func_accuracy.txt");
near_prec<-read.table("results/neighbour_precision_func_accuracy.txt");
ace_prec<-read.table("results/ace_precision_func_accuracy.txt");

#merge the datasets together into a single dataframe
data<-merge(pfam,random_prec,by=0);
data<-merge(data,near_prec,by.x="Row.names",by.y='row.names');
data<-merge(data,ace_prec,by.x="Row.names",by.y='row.names');
data<-merge(data,pfam_meta,by.x="Row.names",by.y='row.names');

#make the variables (columns) of data locally accessible (e.g. 'coverage' instead of 'data$coverage')
attach(data)

#plot coverage vs accuracy
png("figures/pfam_accuracy_NN_vs_coverage.png")
plot(coverage,neighbour_precision,xlab="PFAM Coverage",ylab="PFAM Accuracy (NN Precision)",main="PFAM Accuracy vs Coverage")
#calculate correlation and put it on figure
mtext(sprintf("R2 = %.3f",cor(coverage,neighbour_precision)))
graphics.off()


#if running manually you can call identify so you can click on points and get names
#identify(coverage,neighbour_precision,labels=Row.names)

#Same figure as above but with transposase PFAMs colored red
png("figures/pfam_accuracy_NN_vs_coverage_transposases.png")
plot(coverage,neighbour_precision,xlab="PFAM Coverage",ylab="PFAM Accuracy (NN Precision)",main="Transposases show lower accuracy",pch="*",col="Grey")

####pull out "Transposase" PFAMs
transposases<-unique(c(grep("transposase",pfam_description,ignore.case=TRUE),grep("transposase",Row.names,ignore.case=TRUE)))
integrases<-grep("integrase",pfam_description,ignore.case=TRUE)
mobile_pfams<-unique(c(transposases,integrases))

phage<-unique(c(grep("phage",pfam_description,ignore.case=TRUE),grep("phage",Row.names,ignore.case=TRUE)))

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
png("figures/pfam_accuracy_ace_vs_coverage.png")
plot(coverage,ace_precision,xlab="PFAM Coverage",ylab="PFAM Accuracy (ACE Precision)",main="PFAM Accuracy vs Coverage")
#calculate correlation and put it on figure
mtext(sprintf("R2 = %.3f",cor(coverage,ace_precision)))
graphics.off()

####plot coverage vs accuracy for random
png("figures/pfam_random_accuracy_vs_coverage.png")
plot(coverage,random_precision,xlab="PFAM Coverage",ylab="PFAM Accuracy (Random Precision)",main="PFAM Accuracy (Random) vs Coverage")
#calculate correlation and put it on figure
mtext(sprintf("R2 = %.3f",cor(coverage,random_precision)))
graphics.off()

#plot PFAM coverage distribution
png("figures/pfam_coverage.png")
plot(density(coverage),xlab="PFAM Coverage (# genomes >1 PFAM)",main="Distribution of PFAM Coverage")
graphics.off()

#plot PFAM accuracy distribution
png("figures/PFAM_accuracy_distribution.png")
plot(density(neighbour_precision),col="blue",main="PFAM Accuracy Distribution",xlab="PFAM Accuracy Across Genomes")
lines(density(ace_precision),col="red")
legend("bottom",c("NN","ACE"),col=c("blue","red"),lty=1)
graphics.off()