# predict_tip_states.r
# Author: Dan Knights
# Note: this is in alpha version, code not unit tested.
#
# For help:
# R --slave --vanilla --args -h < predict_tip_states.r
# 
# usage:
# R --slave --vanilla --args -t tree.txt -s tipstates.txt -m KNN -v < predict_tip_states.r
# Neighbors-within-distance method is currently disabled for debugging

library('optparse',warn=FALSE,quiet=TRUE,verbose=FALSE)

# load args
# make option list and parse command line
option_list <- list(
    make_option(c("-t",'--tree_file'), type="character",
        help="tree file [required]."),
    make_option(c("-s", "--states_file"), type="character",
        help="tip states file [required]."),
    make_option(c("-m", "--method"), type="character", default='ACE',
        # help="method: knn (k nearest neighbors), nwd (all neighbors within a fixed branch length distance), or ace (ancestral state reconstruction) [default %default]."),
        help="method: knn (k nearest neighbors) or ace (ancestral state reconstruction) [default %default]."),
    make_option(c("-r", "--ratemodel"), type="character", default='SYM',
        help="rate model for ACE: SYM (symmetrical), ARD (all rates different) [default %default]"),
    make_option(c("--startcol"), type="integer", default=1,
        help="index of start column of tip states file to predict [default predict all]"),
    make_option(c("--endcol"), type="integer", default=NULL,
        help="index of end column of tip states file to predict [default predict all]"),
    make_option(c("-k","--knn_k"), type="integer", default=1,
        help="number of nearest neighbors for KNN [default %default]"),
    # make_option(c("--maxdist"), type="numeric", default=0.01,
    #     help="maximum distance for NWD method [default %default]"),
    make_option(c("--fixed_rate"), type='numeric', default=NULL,
        help="fixed evolutionary rate to prevent ACE from inferring rate matrix [default %default]"),
    make_option(c("-o", "--outdir"), type="character", default='.',
        help="Output directory [default %default]"),
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
        help="print warnings and other additional information")
)
opts <- parse_args(OptionParser(option_list=option_list), args=commandArgs(trailingOnly=TRUE))
dir.create(opts$outdir,showWarnings=FALSE)
library(ape, warn=F, quiet=T, verbose=F)
source('utils.r')
source('tree_utils.r')
source('ace_wrapper.r')
source('prediction.r')

# load data
inputs <- load_tree_and_tipstates(opts$tree_file, opts$states_file)
tree <- inputs$tree
tipstates <- inputs$tipstates
tiptipdists <- inputs$tiptipdists

# if end column not specified, choose last column
if(is.null(opts$endcol)) opts$endcol <- ncol(tipstates)

# yhat is a matrix of predictions, num organisms x num genes
yhat <- matrix(NA, nrow=nrow(tipstates), ncol=(opts$endcol-opts$startcol+1))
rownames(yhat) <- rownames(tipstates)
colnames(yhat) <- colnames(tipstates)[opts$startcol:opts$endcol]

if(opts$verbose){
    cat(sprintf('Tree has %d tips\n',length(tree$tip.label)))
    cat(sprintf('Tipstates has %d rows and %d cols\n',nrow(tipstates), ncol(tipstates)))
    cat(sprintf('Running columns %d-%d.\n', opts$startcol, opts$endcol))
    cat('Method is', opts$method,'\n')
}

# run inference
basefname <- sprintf('method_%s', opts$method)
if(opts$startcol > 1 || opts$endcol < ncol(tipstates))
    basefname <- sprintf('columns_%05d-%05d_%s', opts$startcol, opts$endcol, basefname)
if(opts$method=='ace'){
    library(Matrix, warn=F, quiet=T, verbose=F)
    
    # add rate model to base filename
    basefname <- sprintf('%s_ratemodel_%s', basefname, opts$ratemodel)
    # add fixed rate to base filename
    if(!is.null(opts$fixed_rate)) basefname <- sprintf('%s_fixedrate_%f', basefname, opts$fixed_rate)

    gainrates <- matrix(NA, nrow=nrow(tipstates), ncol=(opts$endcol-opts$startcol+1))
    rownames(gainrates) <- rownames(tipstates)
    colnames(gainrates) <- colnames(tipstates)[opts$startcol:opts$endcol]
    lossrates <- gainrates
    
    cat('Running ACE...\n')
    for(column in opts$startcol:opts$endcol){
        looresults <- loo.ace(tipstates[,column], tree, verbose=opts$verbose, model=opts$ratemodel, fixed.rate=opts$fixed_rate)    
        yhat[,column-opts$startcol + 1] <- looresults$yhat
        if(opts$ratemodel=='SYM'){
            gainrates[,column] <- looresults$rates
        } else {
            gainrates[,column] <- looresults$rates[,1]
            lossrates[,column] <- looresults$rates[,2]
        }
    }
} else if(opts$method=='knn'){
    # add K to base filename
    basefname <- sprintf('%s_K_%d', basefname, opts$knn_k)
    
    cat(sprintf('Running KNN with k=%d...\n', opts$knn_k))
    for(column in opts$startcol:opts$endcol)
        yhat[,column-opts$startcol + 1] <- get.knn.predictions(tipstates[,column], tree, tiptipdists, k=opts$knn_k)
# } else if(opts$method=='nwd'){
#     # add maxdist to base filename
#     basefname <- sprintf('%s_maxdist_%f', basefname, opts$maxdist)
#     cat(sprintf('Running NWD with maxdist=%f...\n',opts$maxdist))
#     
#     for(column in opts$startcol:opts$endcol)
#         yhat[,column-opts$startcol + 1] <- get.knn.predictions.by.distance(tipstates[,column], tree, tiptipdists, distance=opts$maxdist)
#         
} else{
    stop('Unknown inference method.')
}


# save results
fname <- sprintf('%s/predicted_tip_states_%s.txt', opts$outdir, basefname)
sink(fname)
options(digits=5)
cat('Organism\t')
write.table(yhat, sep='\t', quote=FALSE)
sink(NULL)

# save rate matrix if used ACE
if(opts$method=='ACE'){
    fname <- sprintf('%s/inferred_gain_rates_%s.txt', opts$outdir, basefname)
    sink(fname)
    options(digits=5)
    cat('Organism\t')
    write.table(gainrates, sep='\t', quote=FALSE)
    sink(NULL)

    if(opts$ratemodel=='ARD'){
        fname <- sprintf('%s/inferred_loss_rates_%s.txt', opts$outdir, basefname)
        sink(fname)
        options(digits=5)
        cat('Organism\t')
        write.table(lossrates, sep='\t', quote=FALSE)
        sink(NULL)
    }
}