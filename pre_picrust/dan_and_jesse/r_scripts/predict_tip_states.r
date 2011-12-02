# usage:
# -t tree file
# -s tip states file
# -v verbose (optional)
# -m method (KNN, NWD, or ACE)
# -o output directory (must exit)
# --ratemodel (SYM, ARD) (symmetrical, All-Rates-Different)
# --startcol index of start column of tip states file to predict, else predicts all
# --endcol index of end column
# --maxdist maximum distance for NWD (neighbors within distance) method
# -k k for knn
# --fixedrate fixed rate to prevent ACE from inferring rate matrix
# R --slave --vanilla --args -t tree.txt -s tipstates.txt -c 1 -m ACE -v < predict_tip_states.r

source('utils.r')

# load args
arglist <- commandArgs(trailingOnly=TRUE)
treefile <- arglist[which(arglist=='-t')+1]
tipfile <- arglist[which(arglist=='-s')+1]
method <- arglist[which(arglist=='-m')+1]

# optional args
if(sum(arglist=='--startcol')>0){
    startcol <- as.numeric(arglist[grep('--startcol', arglist)+1])
} else startcol <- 1

if(sum(arglist=='--endcol')>0){
    endcol <- as.numeric(arglist[grep('--endcol', arglist)+1])
} else endcol <- NULL

if(sum(arglist=='--maxdist')>0) {
    maxdist <- as.numeric(arglist[which(arglist=='--maxdist')+1])
} else maxdist <- NULL

if(sum(arglist=='--fixedrate')>0){
    fixedrate <- as.numeric(arglist[grep('--fixedrate', arglist)+1])
} else fixedrate <- NULL

if(sum(arglist=='--ratemodel')>0){
    ratemodel <- arglist[grep('--ratemodel', arglist)+1]
} else ratemodel <- 'SYM'

if(sum(arglist=='-k')>0){
    K <- as.numeric(arglist[which(arglist=='-k')+1])
} else K <- NULL

if(sum(arglist=='-o')>0){
    outdir <- arglist[which(arglist=='-o')+1]
} else outdir <- '.'
dir.create(outdir)

verbose <- sum(arglist=='-v') > 0

# load data
inputs <- load_tree_and_tipstates(treefile, tipfile)
tree <- inputs$tree
tipstates <- inputs$tipstates
tiptipdists <- inputs$tiptipdists

# if end column not specified, choose last column
if(is.null(endcol)) endcol <- ncol(tipstates)

# yhat is a matrix of prediction, num organisms x num columns
yhat <- matrix(NA, nrow=nrow(tipstates), ncol=(endcol-startcol+1))
rownames(yhat) <- rownames(tipstates)
colnames(yhat) <- colnames(tipstates)[startcol:endcol]

if(verbose){
    cat(sprintf('Tree has %d tips\n',length(tree$tip.label)))
    cat(sprintf('Tipstates has %d rows and %d cols\n',nrow(tipstates), ncol(tipstates)))
    cat(sprintf('Running columns %d-%d.\n', startcol, endcol))
    cat('Method is', method,'\n')
}



# run inference
basefname <- sprintf('columns_%05d-%05d_method_%s', startcol, endcol, method)
if(method=='ACE'){
    # add fixed rate to base filename
    if(!is.null(fixedrate)) basefname <- sprintf('%s_fixedrate_%f', basefname, fixedrate)
    # add rate model to base filename
    basefname <- sprintf('%s_ratemodel_%s', basefname, ratemodel)

    gainrates <- matrix(NA, nrow=nrow(tipstates), ncol=(endcol-startcol+1))
    rownames(gainrates) <- rownames(tipstates)
    colnames(gainrates) <- colnames(tipstates)[startcol:endcol]
    lossrates <- gainrates
    
    cat('About to run ACE...\n')
    for(column in startcol:endcol){
        looresults <- loo.ace(tipstates[,column], tree, verbose=verbose, model=ratemodel, fixed.rate=fixedrate)    
        yhat[,column-startcol + 1] <- looresults$yhat
        if(ratemodel=='SYM'){
            gainrates[,column] <- looresults$rates
        } else {
#            cat('ARD: ')
#            print(class(looresults$rates))
            gainrates[,column] <- looresults$rates[,1]
            lossrates[,column] <- looresults$rates[,2]
        }
    }
} else if(method=='KNN'){
    # add K to base filename
    basefname <- sprintf('%s_K_%d', basefname, K)
    
    cat('About to run KNN...\n')
    if(is.null(K)) stop("Please specify K.")
    for(column in startcol:endcol)
        yhat[,column-startcol + 1] <- get.knn.predictions(tipstates[,column], tree, tiptipdists, k=K)
} else if(method=='NWD'){
    # add maxdist to base filename
    basefname <- sprintf('%s_maxdist_%f', basefname, maxdist)

    cat('About to run NWD...\n')
    if(is.null(maxdist)) stop("Please specify maxdist.")
    
    for(column in startcol:endcol)
        yhat[,column-startcol + 1] <- get.knn.predictions.by.distance(tipstates[,column], tree, tiptipdists, distance=maxdist)
        
} else{
    stop('Unknown inference method.')
}

# 

# save results
fname <- sprintf('%s/predicted_tip_states_%s.txt', outdir, basefname)
sink(fname)
options(digits=5)
cat('Organism\t')
write.table(yhat, sep='\t', quote=FALSE)
sink(NULL)

# save rate matrix if used ACE
if(method=='ACE'){
    fname <- sprintf('%s/inferred_gain_rates_%s.txt', outdir, basefname)
    sink(fname)
    options(digits=5)
    cat('Organism\t')
    write.table(gainrates, sep='\t', quote=FALSE)
    sink(NULL)

    if(ratemodel=='ARD'){
        fname <- sprintf('%s/inferred_loss_rates_%s.txt', outdir, basefname)
        sink(fname)
        options(digits=5)
        cat('Organism\t')
        write.table(lossrates, sep='\t', quote=FALSE)
        sink(NULL)
    }
}


# how we plotted the "by distance" results
#~ cameans <- apply(cas.by.dist$cas, 2, mean)
#~ nanmeans <- apply(cas.by.dist$nan.rate, 2, mean)
#~ totalmeans <- (1 - (nanmeans))*cameans
#~ plot(cas.by.dist$dists.tried, 1-apply(cas.by.dist$cas, 2, mean), xlab="Distance used", ylab="Classification error",type='b', ylim=ylim)
#~ points(cas.by.dist$dists.tried, apply(cas.by.dist$nan.rate, 2, mean), type='b',pch=15)
#~ points(cas.by.dist$dists.tried, 1-totalmeans, type='b',pch=15,col='blue')
