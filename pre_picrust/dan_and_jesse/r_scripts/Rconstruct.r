library(Matrix)
library(ape)


# note: only uses symmetrical rates
"loo.ace" <- function(states, tree, verbose=FALSE, ...){
    # predicted values
    yhat <- rep(-1,length(states))
    for(tipix in 1:length(states)){
        tipname <- tree$tip.label[tipix]
#~         cat(sprintf('Getting loo estimates for tip %s\n', tipname))
        if(verbose && tipix %% 10 == 0) cat('.')
        # get parent and grandparent's name and distance
        parent.name <- get.parent.name(tree, tipname)
        parent.dist <- get.parent.dist(tree, tipname)
        gparent.name <- get.parent.name(tree, parent.name)
        gparent.dist <- get.parent.dist(tree, parent.name)
        # if this tip was a child of root, use root as the grandparent
        if(is.null(gparent.name)){
            gparent.name <- parent.name
            gparent.dist <- parent.dist
        }
        
        loo.tree <- drop.tip(tree, tipname)
        loo.states <- states[-tipix]

        # if all remaining states are 0 or 1, just predict that 
        if(mean(loo.states)==0){
            yhat[tipix] <- 0
        } else if(mean(loo.states)==1){
            yhat[tipix] <- 1
        } else {
            res <- ace(loo.states, loo.tree, ...)

            # assume parent not in loo tree FIX THIS!!! 
            
            # get rate matrix
            Q <- matrix(res$rates,nrow=2,ncol=2)
            diag(Q) <- -res$rates
            
            # get p = normalized likelihood of PRESENT at gparent
#~             ancestor.name <- get.nearest.ancestor.in.subtree(tree,loo.tree,tipix)
#~             ancestor.ix <- which(loo.tree$node.label==ancestor.name)
            gparent.node.ix <- which(tree$node.label==gparent.name)
            p.gparent <- res$lik.anc[gparent.node.ix,2]
#~             p.ancestor <- res$lik.anc[ancestor.ix,2]
#~             ancestor.dist <- get.ancestor.dist(tree, tipname, ancestor.name)
#~             print(tipname)
#~             print(ancestor.name)
#~             print(ancestor.dist)
            # exp(qt)
            pmatrix <- expm(Q * gparent.dist)
#~             pmatrix <- expm(Q * ancestor.dist)
            
            # prob parent is expectation over PRESENT in pmatrix
            p <- p.gparent * pmatrix[2,2] + (1-p.gparent) * pmatrix[1,2]
#~             p <- p.ancestor * pmatrix[2,2] + (1-p.ancestor) * pmatrix[1,2]
            yhat[tipix] <- p
        }
    }
    if(verbose) cat('\n')
    return(yhat)
}

"get.random.subtree" <- function(tree, fraction){
    ntips <- length(tree$tip.label)
    ndrop <- round((1-fraction) * ntips)
    tipixs <- sample(ntips, ndrop)
    return(drop.tip(tree, tipixs))
}

"get.nearest.ancestor.in.subtree" <- function(tree, subtree, tipix){
    # start with the name of the parent
    tipname <- tree$tip.label[tipix]
    ancestor <- get.parent.name(tree,tipname)
    subtreeix <- which(subtree$node.label==ancestor)
    # while the ancestor isn't in the subtree
    # try the ancestor's parent
    while(length(subtreeix)==0){
        ancestor <- get.parent.name(tree,ancestor)
        if(is.null(ancestor)) return(NULL)
        subtreeix <- which(subtree$node.label==ancestor)
    }
    return(subtree$node.label[subtreeix])
}

"get.parent.name" <- function(tree, nodename){
    # case 1: node is a tip
    ntips <- length(tree$tip.label)
    root.ix <- 1 + ntips
    if(sum(tree$tip.label==nodename) > 0){
        node.ix <- which(tree$tip.label==nodename)
        # root is always the first internal node, return NULL for parent of root
        if(node.ix == root.ix) return(NULL)
        edge.ix <- which(tree$edge[,2]==node.ix)
        parent.ix <- tree$edge[edge.ix,1]
    } else {
        # case 2: node is not a tip
        node.ix <- which(tree$node.label==nodename) + ntips
        # root is always the first internal node, return NULL for parent of root
        if(node.ix == root.ix) return(NULL)
        edge.ix <- which(tree$edge[,2]==node.ix)
        parent.ix <- tree$edge[edge.ix,1]
    }
    parent.name <- tree$node.label[parent.ix - ntips]
    return(parent.name)
}

"get.parent.dist" <- function(tree, nodename){
    # case 1: node is a tip
    ntips <- length(tree$tip.label)
    root.ix <- 1 + ntips
    if(sum(tree$tip.label==nodename) > 0){
        node.ix <- which(tree$tip.label==nodename)
    } else {
        # case 2: node is not a tip
       node.ix <- which(tree$node.label==nodename) + ntips
    }
    if(node.ix == root.ix) return(NULL)
    edge.ix <- which(tree$edge[,2]==node.ix)
    return(tree$edge.length[edge.ix])
}

# note: this uses dist.nodes, which is probably slow!
"get.ancestor.dist" <- function(tree, nodename, ancestorname){
    # case 1: node is a tip
    ntips <- length(tree$tip.label)
    root.ix <- 1 + ntips
    if(sum(tree$tip.label==nodename) > 0){
        node.ix <- which(tree$tip.label==nodename)
    } else {
        # case 2: node is not a tip
       node.ix <- which(tree$node.label==nodename) + ntips
    }
    ancestor.ix <- which(tree$node.label==ancestorname) + ntips
    return(dist.nodes(tree)[node.ix, ancestor.ix])
}


"get.nearest.tip.ix" <- function(tree, tipix){
    return(sort(dists[tipix,],index=T)$ix[-1][1])
}

"get.nearest.tips.ix" <- function(tree, tipix, k=1){
    return(sort(dists[tipix,],index=T)$ix[-1][1:k])
}

"get.nearest.tips.ix.by.distance" <- function(tree, tipix, distance=.03){
    ixs <- which(dists[tipix,] <= distance)
    return(ixs[ixs != tipix])
}

"get.nn.predictions" <- function(states, tree){
    yhat <- rep(-1,length(states))
    for(tipix in 1:length(states)){
        yhat[tipix] <- states[get.nearest.tip.ix(tree, tipix)]
    }
    return(yhat)
}

"get.knn.predictions" <- function(states, tree, k=1){
    yhat <- rep(-1,length(states))
    for(tipix in 1:length(states)){
        yhat[tipix] <- round(mean(states[get.nearest.tips.ix(tree, tipix, k=k)]))
    }
    return(yhat)
}

"get.knn.predictions.by.distance" <- function(states, tree, distance=.03){
    yhat <- rep(-1,length(states))
    for(tipix in 1:length(states)){
        neighbor.ixs <- get.nearest.tips.ix.by.distance(tree, tipix, distance=distance)
        yhat[tipix] <- round(mean(states[neighbor.ixs]))
    }
    return(yhat)
}

# returns matrix of species x gene predictions, same size as tipstates
"predict.all" <- function(tree, tipstates, type=c('ace','knn')[1], verbose=FALSE){
    cas <- rep(0,ncol(tipstates))
    yhat <- matrix(-1,nrow(tipstates), ncol(tipstates))
    for(i in 1:ncol(tipstates)){
        if(verbose) cat(sprintf('Getting L.O.O. predictions for column %d...\n',i))
        if(type=='ace'){
            yhat.i <- loo.ace(tipstates[,i], tree, verbose=verbose, type="discrete", model="SYM", method="ML")
        } else {
            yhat.i <- get.nn.predictions(tipstates[,i], tree)
        }
        yhat[,i] <- as.numeric(yhat.i > .5)
        ca <- sum(yhat[,i] == tipstates[,i])/nrow(tipstates)
        if(verbose) cat(sprintf('...Classification accuracy = %.3f\n', ca))
    }
    return(yhat)
}


"test.different.values.of.k" <- function(tree, tipstates, nKs=100, ngenes=100){
    cas <-matrix(-1,ngenes,nKs); 
    for(geneix in 1:ngenes){
        for(k in 1:nKs){
            y <- tipstates[,geneix]
            yhat <- get.knn.predictions(tipstates[,geneix], tree, k=k)
            ca.gene.k <- sum(yhat == y)/length(yhat); 
            cas[geneix,k] <- ca.gene.k
        }
        if(geneix %% 10) cat(round(10 * geneix/ngenes));
    }
    cat('\n')
    return(cas)
}



"test.different.distances" <- function(tree, tipstates, ndists=100, mindist=0.01, maxdist=max(dists), ngenes=100){
    dists.to.try <- seq(mindist, maxdist,length=ndists)
    cas <-matrix(-1,ngenes,ndists); 
    nan.rate <-matrix(-1,ngenes,ndists); 
    for(geneix in 1:ngenes){
        for(dist.ix in 1:ndists){
            d <- dists.to.try[dist.ix]
            y <- tipstates[,geneix]
            yhat <- get.knn.predictions.by.distance(tipstates[,geneix], tree, distance=d)
            pct.nans <- sum(is.nan(yhat))/length(y)
            nan.ixs <- is.nan(yhat)
            pct.err.within.non.nans <- sum(yhat[!nan.ixs] != y[!nan.ixs])/length(y[!nan.ixs])
            cas[geneix,dist.ix] <- 1-pct.err.within.non.nans
            nan.rate[geneix,dist.ix] <- pct.nans
        }
        if(geneix %% 10) cat(round(10 * geneix/ngenes));
    }
    cat('\n')
    return(list(cas=cas,nan.rate=nan.rate,dists.tried=dists.to.try))
}

# LOAD DATA
#~ if(TRUE){
if(FALSE){
    tree <- read.tree('../validation/ko/data/16S_seqs_nuc_vs_GreenGenesCore_May_5.fasta_aligned_pfiltered_colons_removed.tree.noparalog')
#~     # test set with only 5 genes
#~     tipstates <- read.table('../validation/ko_small_subset/data/kegg_genome_data_April_15_2010_no_ambig_test.fna_vector.txt', row.names=1, header=F,sep='\t')
    # all 4762 genes
    tipstates <- read.table('../validation/ko/data/kegg_genome_data_April_15_2010.fna_vector.txt', row.names=1, header=F,sep='\t')
    for(label in tree$tip.label){if(sum(rownames(tipstates)==label)==0) tree = drop.tip(tree,label)}
    tipstates <- tipstates[tree$tip.label,]
    tree <- multi2di(tree)
    for(i in 1:length(tree$node.label)){if(tree$node.label[i]==""){tree$node.label[i] <- sprintf('inode%d',i)}}
    # fix 0-length edges in tree
    tree$edge.length[tree$edge.length==0] <- 1e-10
    dists <- cophenetic.phylo(x=tree)
}


# how we plotted the "by distance" results
#~ cameans <- apply(cas.by.dist$cas, 2, mean)
#~ nanmeans <- apply(cas.by.dist$nan.rate, 2, mean)
#~ totalmeans <- (1 - (nanmeans))*cameans
#~ plot(cas.by.dist$dists.tried, 1-apply(cas.by.dist$cas, 2, mean), xlab="Distance used", ylab="Classification error",type='b', ylim=ylim)
#~ points(cas.by.dist$dists.tried, apply(cas.by.dist$nan.rate, 2, mean), type='b',pch=15)
#~ points(cas.by.dist$dists.tried, 1-totalmeans, type='b',pch=15,col='blue')

# reads in a tree and tip state matrix, keeps only species in both
"load_tree_and_tipstates" <- function(treefile, tipfile) {
    tree <- read.tree(treefile)
    tipstates <- read.table(tipfile, row.names=1, header=F,sep='\t')
    for(label in tree$tip.label){
        if(sum(rownames(tipstates)==label)==0) 
            tree = drop.tip(tree,label)
    }
    tipstates <- tipstates[tree$tip.label,]
    # convert tree to bifurcating
    tree <- multi2di(tree)
    # add names to new internal nodes
    for(i in 1:length(tree$node.label)){
        if(tree$node.label[i]=="")
            tree$node.label[i] <- sprintf('inode%d',i)
    }
    # fix 0-length edges in tree
    tree$edge.length[tree$edge.length==0] <- 1e-10
    dists <- cophenetic.phylo(x=tree)
    return(list(tree=tree, tipstates=tipstates, tiptipdists=dists))
}
