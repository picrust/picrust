# Predict each tip's state (0 or 1) using ancestral state reconstruction.
# Predictions are "leave-one-out"
# model and method are passed to "ace" function
"loo.ace" <- function(states, tree, verbose=FALSE, 
        type="discrete", model="SYM", method="ML", fixed.rate=NULL){
            
    # predicted values
    yhat <- rep(NA,length(states))
    if(model=='SYM'){
        rates <- rep(NA, length(states))
    } else {
        rates <- matrix(NA, length(states), 2)
    }

    # Loop over tips, hold out each one and predict
    for(tipix in 1:length(states)){
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
        
        # Drop this tip from the tree to get a "leave-one-out" tree
        loo.tree <- drop.tip(tree, tipname)
        loo.states <- states[-tipix]

        # if all remaining states are 0 or 1, just predict that 
        if(mean(loo.states)==0){
            yhat[tipix] <- 0
        } else if(mean(loo.states)==1){
            yhat[tipix] <- 1
        } else {
            if(!is.null(fixed.rate)){
                res <- ace.fixed.rate(loo.states, loo.tree,
                            type=type, model=model, method=method, ip=fixed.rate)
            } else {
                res <- ace(loo.states, loo.tree,
                            type=type, model=model, method=method)
            }

            # create the rate matrix
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
            # (because the parent is dropped from the tree)
            gparent.node.ix <- which(tree$node.label==gparent.name)
            p.gparent <- res$lik.anc[gparent.node.ix,2]
            pmatrix <- expm(Q * gparent.dist)

            # prob parent is expectation over PRESENT in pmatrix
            p <- p.gparent * pmatrix[2,2] + (1-p.gparent) * pmatrix[1,2]
            yhat[tipix] <- p
        }
    }
    if(verbose) cat('\n')
    names(yhat) <- rownames(tipstates)
    return(list(yhat=yhat, rates=rates))
}


# Predict each tip's state (0 or 1) using the mean of the k nearest neighbors
# distances is the output of dist.nodes(tree)
"get.knn.predictions" <- function(states, tree, distances, k=1){
    yhat <- rep(-1,length(states))
    for(tipix in 1:length(states)){
        yhat[tipix] <- mean(states[get.nearest.tips.ix(tree, tipix, distances, k=k)])
    }
    return(yhat)
}

# Predict each tip's state (0 or 1) using the 
# mean of neighbors within a fixed distance
# distances is the output of dist.nodes(tree)
"get.knn.predictions.by.distance" <- function(states, tree, distances, distance=.03){
    yhat <- rep(NA,length(states))
    for(tipix in 1:length(states)){
        neighbor.ixs <- get.nearest.tips.ix.by.distance(tree, tipix, distances, distance=distance)
        if(length(neighbor.ixs)>0)
            yhat[tipix] <- mean(states[neighbor.ixs])
    }
    return(yhat)
}