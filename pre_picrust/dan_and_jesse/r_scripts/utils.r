library(Matrix)
library(ape)


# note: if not using symmetrical rates, returns 0->1 in first column, 1->0 in second
"loo.ace" <- function(states, tree, verbose=FALSE, 
        type="discrete", model="SYM", method="ML", fixed.rate=NULL){
    # predicted values
    yhat <- rep(NA,length(states))
    if(model=='SYM'){
        rates <- rep(NA, length(states))
    } else {
        rates <- matrix(NA, length(states), 2)
    }

    for(tipix in 1:2){#:length(states)){
        tipname <- tree$tip.label[tipix]
        if(verbose) cat(sprintf('Getting loo estimates for tip %s\n', tipname))
        if(verbose && tipix %% 10 == 0) cat(tipix, ' ')
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
        
            if(!is.null(fixed.rate))
                res <- ace.fixed.rate(loo.states, loo.tree,
                            type=type, model=model, method=method, ip=fixed.rate)
            else
                res <- ace(loo.states, loo.tree,
                            type=type, model=model, method=method)

            # assume parent not in loo tree FIX THIS!!! 
            
            # get rate matrix
            if(model=='SYM'){
                Q <- matrix(res$rates, nrow=2, ncol=2)
                diag(Q) <- -res$rates
                rates[tipix] <- Q[1,2]
            } else {
                Q <- matrix(0, 2, 2)
                Q[2,1] <- res$rates[1]
                Q[1,2] <- res$rates[2]
                Q[1,1] <- -Q[2,1]
                Q[2,2] <- -Q[1,2]
                # 0->1 rates in first column
                rates[tipix,1] <- Q[1,2]
                rates[tipix,2] <- Q[2,1]
            }            
            
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
            
# cat("Prob. grandparent:", p.gparent, '\n')
# cat("Rate matrix: "); print(Q)
# cat("Branch length: "); print(gparent.dist)
# cat("pmatrix: "); print(pmatrix)



            # prob parent is expectation over PRESENT in pmatrix
            p <- p.gparent * pmatrix[2,2] + (1-p.gparent) * pmatrix[1,2]
#~             p <- p.ancestor * pmatrix[2,2] + (1-p.ancestor) * pmatrix[1,2]
            yhat[tipix] <- p
            # print(yhat[tipix])
        }
    }
    if(verbose) cat('\n')

    names(yhat) <- rownames(tipstates)

    return(list(yhat=yhat, rates=rates))
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


"get.nearest.tip.ix" <- function(tree, tipix, distances){
    return(sort(distances[tipix,],index=T)$ix[-1][1])
}

"get.nearest.tips.ix" <- function(tree, tipix, distances, k=1){
    return(sort(distances[tipix,],index=T)$ix[-1][1:k])
}

"get.nearest.tips.ix.by.distance" <- function(tree, tipix, distances, distance=.03){
    ixs <- which(distances[tipix,] <= distance)
    return(ixs[ixs != tipix])
}

"get.nn.predictions" <- function(states, tree, distances){
    yhat <- rep(-1,length(states))
    for(tipix in 1:length(states)){
        yhat[tipix] <- states[get.nearest.tip.ix(tree, tipix, distances)]
    }
    return(yhat)
}

"get.knn.predictions" <- function(states, tree, distances, k=1){
    
    yhat <- rep(-1,length(states))
    for(tipix in 1:length(states)){
        yhat[tipix] <- mean(states[get.nearest.tips.ix(tree, tipix, distances, k=k)])
    }
    return(yhat)
}

"get.knn.predictions.by.distance" <- function(states, tree, distances, distance=.03){
    yhat <- rep(NA,length(states))
    for(tipix in 1:length(states)){
        neighbor.ixs <- get.nearest.tips.ix.by.distance(tree, tipix, distances, distance=distance)
        if(length(neighbor.ixs)>0)
            yhat[tipix] <- mean(states[neighbor.ixs])
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
            yhat.i <- get.knn.predictions(tipstates[,i], tree, k=10)
        }
        yhat[,i] <- yhat.i
        ca <- sum(as.numeric(yhat[,i] > .5) == tipstates[,i])/nrow(tipstates)
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


# reads in a tree and tip state matrix, keeps only species in both
"load_tree_and_tipstates" <- function(treefile, tipfile, min.branch.length=1e-10) {
    tree <- read.tree(treefile)
    tipstates <- read.table(tipfile, row.names=1, header=T)
    # Drop tips that aren't in the tipstates table
    for(label in tree$tip.label){
        if(sum(rownames(tipstates)==label)==0) 
            tree = drop.tip(tree,label)
    }
    
    # drop data rows that aren't in tree
    tipstates <- tipstates[tree$tip.label,]
    
    # convert tree to bifurcating
    tree <- multi2di(tree)
    # add names to new internal nodes
    for(i in 1:length(tree$node.label)){
        if(tree$node.label[i]=="")
            tree$node.label[i] <- sprintf('inode%d',i)
    }
    # fix 0-length edges in tree
    tree$edge.length[tree$edge.length==0] <- min.branch.length
    dists <- cophenetic.phylo(x=tree)
    return(list(tree=tree, tipstates=tipstates, tiptipdists=dists))
}


ace.fixed.rate <- function(x, phy, type = "continuous", method = "ML", CI = TRUE,
                model = if (type == "continuous") "BM" else "ER",
                scaled = TRUE, kappa = 1, corStruct = NULL, ip = 0.1)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo".')
    if (is.null(phy$edge.length))
        stop("tree has no branch lengths")
    type <- match.arg(type, c("continuous", "discrete"))
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if (nb.node != nb.tip - 1)
      stop('"phy" is not rooted AND fully dichotomous.')
    if (length(x) != nb.tip)
      stop("length of phenotypic and of phylogenetic data do not match.")
    if (!is.null(names(x))) {
        if(all(names(x) %in% phy$tip.label))
          x <- x[phy$tip.label]
        else warning("the names of 'x' and the tip labels of the tree do not match: the former were ignored in the analysis.")
    }
    obj <- list()
    if (kappa != 1) phy$edge.length <- phy$edge.length^kappa
    if (type == "continuous") {
        if (method == "pic") {
            if (model != "BM")
              stop('the "pic" method can be used only with model = "BM".')
            ## See pic.R for some annotations.
            phy <- reorder(phy, "pruningwise")
            phenotype <- numeric(nb.tip + nb.node)
            phenotype[1:nb.tip] <- if (is.null(names(x))) x else x[phy$tip.label]
            contr <- var.con <- numeric(nb.node)
            ans <- .C("pic", as.integer(nb.tip), as.integer(nb.node),
                      as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]),
                      as.double(phy$edge.length), as.double(phenotype),
                      as.double(contr), as.double(var.con),
                      as.integer(CI), as.integer(scaled),
                      PACKAGE = "ape")
            obj$ace <- ans[[6]][-(1:nb.tip)]
            names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
            if (CI) {
                se <- sqrt(ans[[8]])
                CI95 <- matrix(NA, nb.node, 2)
                CI95[, 1] <- obj$ace + se * qnorm(0.025)
                CI95[, 2] <- obj$ace - se * qnorm(0.025)
                obj$CI95 <- CI95
            }
        }
        if (method == "ML") {
            if (model == "BM") {
                tip <- phy$edge[, 2] <= nb.tip
                dev.BM <- function(p) {
                    if (p[1] < 0) return(1e100) # in case sigma² is negative
                    x1 <- p[-1][phy$edge[, 1] - nb.tip]
                    x2 <- numeric(length(x1))
                    x2[tip] <- x[phy$edge[tip, 2]]
                    x2[!tip] <- p[-1][phy$edge[!tip, 2] - nb.tip]
                    -2 * (-sum((x1 - x2)^2/phy$edge.length)/(2*p[1]) -
                          nb.node * log(p[1]))
                }
                out <- nlm(function(p) dev.BM(p),
                           p = c(1, rep(mean(x), nb.node)), hessian = TRUE)
                obj$loglik <- -out$minimum / 2
                obj$ace <- out$estimate[-1]
                names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
                se <- sqrt(diag(solve(out$hessian)))
                obj$sigma2 <- c(out$estimate[1], se[1])
                se <- se[-1]
                if (CI) {
                    CI95 <- matrix(NA, nb.node, 2)
                    CI95[, 1] <- obj$ace + se * qt(0.025, nb.node)
                    CI95[, 2] <- obj$ace - se * qt(0.025, nb.node)
                    obj$CI95 <- CI95
                }
            }
        }
        if (method == "GLS") {
            if (is.null(corStruct))
              stop('you must give a correlation structure if method = "GLS".')
            if (class(corStruct)[1] == "corMartins")
              M <- corStruct[1] * dist.nodes(phy)
            if (class(corStruct)[1] == "corGrafen")
              phy <- compute.brlen(attr(corStruct, "tree"),
                                   method = "Grafen",
                                   power = exp(corStruct[1]))
            if (class(corStruct)[1] %in% c("corBrownian", "corGrafen")) {
                dis <- dist.nodes(attr(corStruct, "tree"))
                MRCA <- mrca(attr(corStruct, "tree"), full = TRUE)
                M <- dis[as.character(nb.tip + 1), MRCA]
                dim(M) <- rep(sqrt(length(M)), 2)
            }
            varAY <- M[-(1:nb.tip), 1:nb.tip]
            varA <- M[-(1:nb.tip), -(1:nb.tip)]
            V <- corMatrix(Initialize(corStruct, data.frame(x)),
                           corr = FALSE)
            invV <- solve(V)
            o <- gls(x ~ 1, data.frame(x), correlation = corStruct)
            GM <- o$coefficients
            obj$ace <- drop(varAY %*% invV %*% (x - GM) + GM)
            names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
            if (CI) {
                CI95 <- matrix(NA, nb.node, 2)
                se <- sqrt((varA - varAY %*% invV %*% t(varAY))[cbind(1:nb.node, 1:nb.node)])
                CI95[, 1] <- obj$ace + se * qnorm(0.025)
                CI95[, 2] <- obj$ace - se * qnorm(0.025)
                obj$CI95 <- CI95
            }
        }
    } else { # type == "discrete"
        if (method != "ML")
          stop("only ML estimation is possible for discrete characters.")
        if (!is.factor(x)) x <- factor(x)
        nl <- nlevels(x)
        lvls <- levels(x)
        x <- as.integer(x)
        if (is.character(model)) {
            rate <- matrix(NA, nl, nl)
            if (model == "ER") np <- rate[] <- 1
            if (model == "ARD") {
                np <- nl*(nl - 1)
                rate[col(rate) != row(rate)] <- 1:np
            }
            if (model == "SYM") {
                np <- nl * (nl - 1)/2
                sel <- col(rate) < row(rate)
                rate[sel] <- 1:np
                rate <- t(rate)
                rate[sel] <- 1:np
            }
        } else {
            if (ncol(model) != nrow(model))
              stop("the matrix given as `model' is not square")
            if (ncol(model) != nl)
              stop("the matrix `model' must have as many rows
as the number of categories in `x'")
            rate <- model
            np <- max(rate)
        }
        index.matrix <- rate
        tmp <- cbind(1:nl, 1:nl)
        index.matrix[tmp] <- NA
        rate[tmp] <- 0
        rate[rate == 0] <- np + 1 # to avoid 0's since we will use this as numeric indexing

        liks <- matrix(0, nb.tip + nb.node, nl)
        TIPS <- 1:nb.tip
        liks[cbind(TIPS, x)] <- 1
        phy <- reorder(phy, "pruningwise")
        
        Q <- matrix(0, nl, nl)
        dev <- function(p, output.liks = FALSE) {
            if (any(is.nan(p)) || any(is.infinite(p))) return(1e50)
            ## from Rich FitzJohn:
            comp <- numeric(nb.tip + nb.node) # Storage...
            Q[] <- c(p, 0)[rate]
            diag(Q) <- -rowSums(Q)
            for (i  in seq(from = 1, by = 2, length.out = nb.node)) {
                j <- i + 1L
                anc <- phy$edge[i, 1]
                des1 <- phy$edge[i, 2]
                des2 <- phy$edge[j, 2]
                v.l <- matexpo(Q * phy$edge.length[i]) %*% liks[des1, ]
                v.r <- matexpo(Q * phy$edge.length[j]) %*% liks[des2, ]
                v <- v.l * v.r
                comp[anc] <- sum(v)
                liks[anc, ] <- v/comp[anc]
            }
            if (output.liks) return(liks[-TIPS, ])
            dev <- -2 * sum(log(comp[-TIPS]))
            if (is.na(dev)) Inf else dev
        }
        # cat('about to run optimization of rates...\n')
        # cat('there are', np, 'rates to estimate.\n')
        # SKIPPING THEIR optimization step
        # out <- nlminb(rep(ip, length.out = np), function(p) dev(p),
        #               lower = rep(0, np), upper = rep(1e50, np))
        # obj$loglik <- -out$objective/2
        # obj$rates <- out$par
        obj$loglik <- 0
        obj$rates <- rep(ip, length.out = np)
        # cat('finished optimization of rates, rate is', ip, '\n')

        SKIPSE=TRUE
        if(!SKIPSE){
            oldwarn <- options("warn")
            options(warn = -1)
            h <- nlm(function(p) dev(p), p = obj$rates, iterlim = 1,
                     stepmax = 0, hessian = TRUE)$hessian
            options(oldwarn)
            if (any(h == 0))
              warning("The likelihood gradient seems flat in at least one dimension (gradient null):\ncannot compute the standard-errors of the transition rates.\n")
            else obj$se <- sqrt(diag(solve(h)))
            obj$index.matrix <- index.matrix
        }
        # cat('about to get CI...\n')
        if (CI) {
            obj$lik.anc <- dev(obj$rates, TRUE)
            colnames(obj$lik.anc) <- lvls
        }
        # cat('just got CI.\n')
    }
    obj$call <- match.call()
    class(obj) <- "ace"
    obj
}

