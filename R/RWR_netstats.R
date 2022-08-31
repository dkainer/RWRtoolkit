
# Given a multiplex object, generate a similarity matrix of pairwise similarity
# of the layers in that multiplex
overlapSimilarityMultiplex <- function(mpo){
  
  n = mpo$Number_of_Layers
  mat <- matrix(0,nrow=n, ncol=n)
  
  for(i in 1:n){
    for(j in i:n){
      cat(names(mpo)[i],names(mpo)[j],"\n")
      g.int   <- graph.intersection(mpo[[i]],mpo[[j]])
      g.union <- graph.union(mpo[[i]],mpo[[j]])
      mat[i,j] <- ecount(g.int)/ecount(g.union)
      mat[j,i] <- mat[i,j]
    }
  }
  rownames(mat) = names(mpo)[1:n]
  colnames(mat) = names(mpo)[1:n]
  return(mat)
  #return(sum(vect)/length(multiplex@layers))
}

# Given a multiplex object AND an external network file, 
# calculate the similarity of each multiplex layer to that external network
overlapSimilarityMultiplexLayer <- function(mpo, nwfile=NULL){
  if(is.null(nwfile))
  {
    message("you need to supply the filename of the network layer to test against the multiplex!")
    break
  }
  nw      <- data.table::fread(nwfile, col.names = c("from","to","weight"), select=1:3)
  nw.g    <- graph_from_data_frame(nw, directed = F)
  
  vect <- vector(length = mpo$Number_of_Layers)
  for(i in 1:mpo$Number_of_Layers){
    print(names(mpo)[i])
    g.int   <- graph.intersection(mpo[[i]],nw.g)
    g.union <- graph.union(mpo[[i]],nw.g)
    vect[i] <- ecount(g.int)/ecount(g.union)
  }
  names(vect) <- names(mpo)
  return(vect)
}

# Given to network files, calculate their similarity
overlapSimilarityLayerLayer <- function(nwfile1=NULL, nwfile2=NULL){
  if(is.null(nwfile1) || is.null(nwfile2))
  {
    message("you need to supply the filenames of the network layers to test against the each other!")
    break
  }
  nw      <- data.table::fread(nwfile1, col.names = c("from","to","weight"), select=1:3)
  nw.g1    <- graph_from_data_frame(nw, directed = F)
  nw      <- data.table::fread(nwfile2, col.names = c("from","to","weight"), select=1:3)
  nw.g2    <- graph_from_data_frame(nw, directed = F)
  
  g.int   <- graph.intersection(nw.g1, nw.g2)
  g.union <- graph.union(nw.g1, nw.g2)
  return( ecount(g.int)/ecount(g.union) )
}


# Given a multiplex object and a network file, calculate a weighted score of the
# overlap between each layer in the multiplex and the network file
overlapScoreMultiplex <- function(mpo, testfile = NULL)
{
  if(is.null(testfile))
  {
    message("you need to supply the filename of the network layer to test against the multiplex!")
    break
  }
  
  test      <- data.table::fread(testfile, col.names = c("from","to","weight"), select=1:3)
  test.g    <- graph_from_data_frame(test, directed = F)
  cat("nodes in test network:", vcount(test.g), "\n")
  cat("edges in test network:", ecount(test.g), "\n")
  
  #  for(i in 1:mpo$Number_of_Layers){
  scores <- foreach(i = 1:mpo$Number_of_Layers, .combine=c) %do% {  
    print(names(mpo)[i])
    cat("nodes in", names(mpo)[i], ":", vcount(mpo[[i]]), "\n")
    cat("edges in", names(mpo)[i], ":", ecount(mpo[[i]]), "\n")
    
    g.int   <- graph.intersection(test.g, mpo[[i]])
    g.union <- graph.union(test.g, mpo[[i]])
    cat("number of intersecting edges:",ecount(g.int),"\n")
    cat("% intersect of test network:", 100*ecount(g.int)/ecount(test.g),"\n")
    cat("% intersect of",names(mpo)[i],":", 100*ecount(g.int)/ecount(mpo[[i]]),"\n")
    cat("jaccard:", ecount(g.int)/ecount(g.union),"\n")
    score <- sum(E(g.int)$weight_1) / ecount(mpo[[i]])
    cat("normalized score per", names(mpo)[i]," edge:", score, "\n-----------\n" )
    score
  }
  
  return(scores)
}

# Given a flist file and a network file, calculate a weighted score of the
# overlap between each layer in the flist and the network file
overlapScoreFlist <- function(flist=NULL, testfile = NULL)
{
  library(foreach)
  library(iterators)
  if(is.null(nwfile) || is.null(flist))
  {
    message("you need to supply the filename of the test network and the filename of the flist!")
    break
  }
  test      <- data.table::fread(testfile, col.names = c("from","to","weight"), select=1:3)
  test.g    <- graph_from_data_frame(test, directed = F)
  cat("nodes in test network:", vcount(test.g), "\n")
  cat("edges in test network:", ecount(test.g), "\n")
  
  flistdf <- data.table::fread(flist, header=F, col.names = c("nwfile","nwname","nwgroup"), select = 1:3)
  
  cat("layers to test against:\n")
  print(flistdf$nwname)
  
  scores <- foreach(d=iter(flistdf, by='row')) %do%
    {
      print(d$nwname)
      nw      <- data.table::fread(d$nwfile, col.names = c("from","to","weight"), select=1:3)
      nw.g    <- graph_from_data_frame(nw, directed = F)
      cat("nodes in", d$nwname, ":", vcount(nw.g), "\n")
      cat("edges in", d$nwname, ":", ecount(nw.g), "\n")
      
      g.int   <- graph.intersection(nw.g, test.g)
      g.union <- graph.union(nw.g, test.g)
      cat("number of intersecting edges:",ecount(g.int),"\n")
      cat("% intersect of test network:", 100*ecount(g.int)/ecount(test.g),"\n")
      cat("% intersect of",d$nwname,":", 100*ecount(g.int)/ecount(nw.g),"\n")
      cat("jaccard:", ecount(g.int)/ecount(g.union),"\n")
      score <- sum(E(g.int)$weight_1) / ecount(test.g)
      cat("normalized score per", d$nwname," edge:", score, "\n-----------\n" )
      score
    }
  return(scores)
}


# Same as overlapScoreMultiplex but it also makes a Tau vector so the user 
# can supply Tau for other RWR tools.
getTau <- function(mpo, testfile=NULL)
{
  if(is.null(testfile))
  {
    message("you need to supply the filename of the network layer to test against the multiplex!")
    break
  }
  scores <- overlapScoreMultiplex(mpo = mpo, testfile = testfile)
  if(!is.null(scores)) {
    return( mpo$Number_of_Layers*(scores/sum(scores)) )
  }
}

# merge down all layers in a multiplex object, but don't aggregate the edges (i.e. keep all edges).
merged_with_all_edges <- function(mpo)
{
  # merge the subnetworks into one big network
  cat("merging",mpo$Number_of_Layers," network layers down...\n")
  nl        <- mpo$Number_of_Layers
  nw.dflist <- lapply(mpo[1:nl], igraph::as_data_frame)
  nw.df     <- dplyr::bind_rows(nw.dflist)
  nw.df     <- nw.df %>% dplyr::group_by(type) %>% dplyr::mutate(weightnorm = weight/sum(weight))
  nw.merged <- graph_from_data_frame(nw.df , directed = FALSE)
  print("merging network layers down...DONE")
  cat("Merged network has", ecount(nw.merged), "edges and", vcount(nw.merged),"rows\n" )
  return(nw.merged)
}

# merge down all layers and aggregate multi-edges into one edge where 
# edge weight is the number of layers in which the two nodes are connected.
merged_with_edgecounts <- function(mpo, inv=FALSE){
  
  n = mpo$Number_of_Nodes_Multiplex
  A <- as( matrix(0, ncol=n, nrow=n), "dgCMatrix")
  for(i in 1:mpo$Number_of_Layers)
  {
    #as_adjacency_matrix
    A <- A + as_adjacency_matrix(mpo[[i]], sparse=TRUE)
  }
  aggr <- igraph::graph_from_adjacency_matrix(A, mode="undirected", weighted=TRUE, diag=FALSE)
  aggr <- set.vertex.attribute(aggr,name="name",index=V(aggr),value=mpo$Pool_of_Nodes)
  if(inv){
    aggr <- set.edge.attribute(aggr,"weight",index = E(aggr),value=(1/get.edge.attribute(aggr,"weight",index = E(aggr))))
  }
  
  return(aggr)
}

# Exclusivity is proportion of all edges in the multiplex that are found in only one layer.
# Here we show how many edges are found in 1 layer, 2 layers, 3 layers etc...
exclusivity <- function(mpo){
  merged <- merged_with_edgecounts(mpo)  # this gives us the number of edges between each node pair
  for(i in 1:mpo$Number_of_Layers)
  {
    exc <- round( sum(E(merged)$weight==i)/ecount(merged),4 )
    cat("proportion of all edges found in ",i," layers: ", exc,"\n")
  }
}

