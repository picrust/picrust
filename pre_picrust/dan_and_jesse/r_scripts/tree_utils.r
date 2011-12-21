# Get the name of the parent node of the given node
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

# Get the distance to the parent node of the given node
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

# Get the distance to the parent node of the given node
"get.ancestor.dist" <- function(tree, nodename, ancestorname){
    # case 1: node is a tip
    ntips <- length(tree$tip.label)
    if(any(tree$tip.label==nodename)){
        node.ix <- which(tree$tip.label==nodename)
    } else {
        # case 2: node is not a tip
       node.ix <- which(tree$node.label==nodename) + ntips
    }
    ancestor.ix <- which(tree$node.label==ancestorname) + ntips
    return(dist.nodes(tree)[node.ix, ancestor.ix])
}

# Get the index of the nearest tip to the given tip
# distances is the output of dist.nodes(tree)
"get.nearest.tip.ix" <- function(tree, tipix, distances){
    return(sort(distances[tipix,],index=T)$ix[-1][1])
}

# Get the index of the k nearest tips to the given tip
# distances is the output of dist.nodes(tree)
"get.nearest.tips.ix" <- function(tree, tipix, distances, k=1){
    return(sort(distances[tipix,],index=T)$ix[-1][1:k])
}

# Get the index of all tips within a fixed distance of the given tip
# distances is the output of dist.nodes(tree)
"get.nearest.tips.ix.by.distance" <- function(tree, tipix, distances, distance=.03){
    ixs <- which(distances[tipix,] <= distance)
    return(ixs[ixs != tipix])
}

