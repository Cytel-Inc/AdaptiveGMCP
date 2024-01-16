
#Generate weights for all intersection hypothesis
genWeights <- function(w,g,HypothesisName)
{
  n <- length(w)
  basket <- startIntersects <- InterSects <- list('Intersect'=rep(1, n))
  startNodes <- Nodes <-list('Node'=w)
  startEdeges <- Edges <-list('Edges'=g)

  for(genIDX in 1:(n-1))
  {
    nextIntersects <- list()
    nextNodes <- list()
    nextEdeges <- list()

    for (idx in 1:length(startIntersects)) {
      newBranch <- genNextBranch(vParent = startIntersects[[idx]],
                                w=startNodes[[idx]],
                                g=startEdeges[[idx]],
                                basket=basket)
      if(!length(newBranch)==0){
        nIntersect <- newBranch[seq(1,length(newBranch),by=3)] #every 1st/3 entry is for intersect
        nNodes <- newBranch[seq(2,length(newBranch),by=3)] #every 2nd/3 entry is for nodes
        nEdges <- newBranch[seq(3,length(newBranch),by=3)] #every 3rd/3 entry is for edeges

        basket <- append(basket, nIntersect)
        nextIntersects <- append(nextIntersects, nIntersect)
        nextNodes <- append(nextNodes, nNodes)
        nextEdeges <- append(nextEdeges, nEdges)
      }
    }
    startIntersects <- nextIntersects
    startNodes <- nextNodes
    startEdeges <- nextEdeges

    InterSects <- append(InterSects, nextIntersects)
    Nodes <- append(Nodes, nextNodes)
    Edges <- append(Edges, nextEdeges)
  }
  WH <- data.frame(cbind(do.call(rbind, InterSects), do.call(rbind, Nodes)))
  row.names(WH) <- NULL
  IntersectIDX <- as.vector(apply(WH[,1:(ncol(WH)/2)], 1, function(x){paste(x,collapse = '')}))

  if(!missing(HypothesisName))
  {
    colnames(WH) <- c(HypothesisName,
                      paste('Weight',get_numeric_part(HypothesisName),sep=''))
    return(
      list('IntersectionWeights'=WH,
           'Edges'=Edges,
           'IntersectIDX'= IntersectIDX,
           'HypothesisName'=HypothesisName))
  }else
  {
    return(
      list('IntersectionWeights'=WH,
         'Edges'=Edges,
         'IntersectIDX'= IntersectIDX
         ))
  }
}



#Generate new branches form a parent node with the updated weights and transition matrices
genNextBranch <- function(vParent,w,g,basket)
{
  J <- which(vParent != 0)
  newBranch <- list()
  for(j in J)
  {
    b <- vParent
    b[j] <- 0
    bIDX <- paste(b,collapse  = '')
    basketIDX <- lapply(1:length(basket),
                        function(x) paste(basket[[x]],collapse  = ''))

    if(!any(basketIDX==bIDX))
    {
      intHypo <- b
      newGraph <- getUpdatedGraph(w = w,g = g,zero_idx = j)
      newBranch <- append(newBranch,
                          list('Intersect'=b,
                               'Node'=newGraph$Node,
                               'Edges'=newGraph$Edges))
    }
  }
  newBranch
}


#generate the new weights from the parent node and edges and generate new nodes for later computations
getUpdatedGraph <- function(w, g, zero_idx)
{
  w1 <- rep(0, length(w))
  for(i in 1:length(1:length(w)))
  {
    w1[i] <- ifelse(is.na(w[i]) || i == zero_idx,
                    0,
                    w[i] + g[zero_idx,i]*w[zero_idx])
  }
  list("Node"=w1,
       "Edges"=UpdateTransgatrix(g, zero_idx))
}



#Update Transition matrix
UpdateTransgatrix <- function(g, zero_idx)
{
  g1 = g
  for(l in 1:nrow(g))
  {
    for (k in 1:ncol(g)) {
      if(l == zero_idx || k == zero_idx || is.na(g[l,k]))
      {
        g1[l,k] <- 0
      }else
      {
        g1[l,k] = ifelse(k!=l & g[l,zero_idx]*g[zero_idx,l] < 1,
                         (g[l,k] + g[l,zero_idx]*g[zero_idx,k])/(1-g[l,zero_idx]*g[zero_idx,l]),
                         0)

      }
    }
  }
  g1
}





