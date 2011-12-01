#!/usr/bin/env Rscript
#library(hexbin)

library(sciplot)

data<-read.table("results/subsystem_neighbour_precision_func_all_accuracy.txt")
names(data)<-c("gp_id","SS","Accuracy","NN")

metadata<-read.delim("data/subsystem_classes.txt")

all_data<-merge(data,metadata,by.x="SS",by.y='row.names')

cat_list<-unique(all_data$category)

acc_max=max(all_data$Accuracy)
bin_size<-0.05
num_bins<-ceiling(max(all_data$NN)/bin_size)

for(i in 1:length(cat_list)){
png(paste("figures/categories/",cat_list[i],"_binned.png",sep=""))

xy<-all_data[all_data$category==cat_list[i],c('NN','Accuracy')]

j<-0
bin_categories<-numeric(0)
while(j < max(xy$NN)){
	bin_categories[xy$NN >=j & xy$NN <j+bin_size]<-j
	j<-j+bin_size
}
bargraph.CI(bin_categories,xy$Accuracy,err.width=0.05,ylim=c(0,acc_max),xlim=c(0,num_bins),axis.lty=1,main=cat_list[i],ylab="Accuracy (NN precision per subsystem)",xlab=paste("NN 16s tree distance (binned every",bin_size,")"))
graphics.off()
}
quit(save="no")

for(i in 1:length(cat_list)){
png(paste("figures/categories/",cat_list[i],"_hexbin.png",sep=""))

bin_ave<-0
bins<-0
k<-1
j<-0
while(j < max(xy$NN)){
	bin_ave[k]<-mean(xy$Accuracy[xy$NN >=j & xy$NN <j+bin_size])
	bins[k]<-j
	k<-k+1
	j<-j+bin_size
}

xy<-all_data[all_data$category==cat_list[i],c('NN','Accuracy')]
plot(hexbin(xy),main=cat_list[i],ylab="Accuracy (NN precision per subsystem)",xlab="NN 16s tree distance")
#mtext(sprintf("R2 = %.3f",cor(xy$NN,xy$Accuracy)))
graphics.off()
}
