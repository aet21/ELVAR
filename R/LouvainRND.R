#' @title
#' Louvain algorithm detects communities from a graph network
#'
#' @aliases LouvainRND
#'
#' @description
#' This function optimizes (maximize) the modularity function for identifying communities from a graph network
#' following Louvain-like method. Louvain implementation in igraph R package gives a unique network partition
#' from a graph network on every run due to sequential selection of node for local moving during optimization.
#' But this Louvain implementation may generate slightly different network partitions from a graph network on
#' every run because we considered random (nonsequential) selection of node for local moving during optimization.
#'
## Use LouvainRND() function to run the code...
#'
## Input parameters......
#' @param G
#' An igraph object
#' @param resolution
#' Tune the size of the communities. Default is 1. Community size decreases with the increment of resolution value
#' @param threshold
#' If change in the quality function (score) is less than the threshold, iteration loop will be terminated
#'
## Output parameters......
#' @return
#' It is a list consists of
#' CommunityMembers:: is a vector where names are the nodes and the values are the communities belong to the nodes, and
#' ModularityVal:: Modularity value
#'
#' @references
#'
#' @examples
#'
#' @import dplyr
#' @import igraph
#'
#' @export
#'
#'
LouvainRND <- function(G, resolution= 1, threshold= 0.0000001){
  ###

  # Check the input graph is a graph object or not
  if (!is_igraph(G)){
    stop("Not a graph object")
  }
  # Check the input graph is a directed graph or not
  if (is.directed(G)){
    stop("It's a directed graph object")
  }
  # Add edge weight if the input graph isn't weighted
  if(!is_weighted(G)){
    E(G)$weight <- 1
  }
  # Extract vertex name if available in the input graph
  vertexN.v <- NULL
  if(!is.null(names(V(G)))){
    vertexN.v <- names(V(G))
    #
    G <- delete_vertex_attr(G, "name") # Delete vertex "name" attribute
  }
  # Run Louvain optimization
  OneLevelSumm <- Lone_levelFirst(G, resolution, threshold)
  #
  mod <- OneLevelSumm$Modu
  CommInfo <- OneLevelSumm$clusterInfo
  Grph_status <- OneLevelSumm$GrpSta
  Grph_status$commu <- CommInfo
  #
  current_grp <- Lgen_graph(G, CommInfo)
  #
  while(TRUE){
    OneLevelSumm <- Lone_levelSecond(current_grp, resolution, threshold, Grph_status)
    #
    new_mod <- OneLevelSumm$Modu
    CommInfo <- OneLevelSumm$clusterInfo
    Grph_status$commu <- CommInfo

    if((new_mod - mod) <= threshold){

      members <- rep(1:length(CommInfo), lengths(CommInfo))
      names(members) <- unlist(CommInfo)
      members <- members[sort(as.numeric(names(members)), index.return= TRUE)$ix]
      if(!is.null(vertexN.v)){
        names(members) <- vertexN.v
        }
      return(list(CommunityMembers= members, ModularityVal= new_mod))
    }
    #
    mod <- new_mod
    current_grp <- Lgen_graph(G, CommInfo)
   }
}

#### Auxiliary functions
#'
#####****$$$$$$$$$$$$***Graph Partition***First**$$$$$$$$$$$***#####
Lone_levelFirst <- function(G, resolution, threshold){
  #
  GrpAdj <- as.matrix(get.adjacency(G, attr = "weight", sparse = T))
  diag(GrpAdj) <- diag(GrpAdj)*2
  out_degree <- rowSums(GrpAdj)
  in_degree <- out_degree
  norm <- 1/sum(out_degree)
  #
  vlist <- V(G)
  ini_partition <- as.list(vlist)
  node2com <- as.vector(vlist)
  names(node2com) <- seq(length(node2com))
  inner_partition <- ini_partition
  #
  Grph_status <- list(Grp = G, ODeg= out_degree, IDeg= in_degree, Norm= norm, commu= ini_partition)
  #
  degrees.v <- out_degree
  print(sum(degrees.v))
  Stot <- degrees.v
  #
  cur_mod <- LNModularity(Grph_status, resolution, norm)
  cur_Grph_status <- Grph_status
  new_mod <- cur_mod
  #
  diag(GrpAdj) <- 0
  # nbrs <- lapply(1:length(node2com), function(x) cbind(as.numeric(which(GrpAdj[x, ] > 0)), GrpAdj[x, GrpAdj[x, ] > 0]))
  nbrs <- lapply(1:length(node2com), function(x) {Selidx <- as.numeric(which(GrpAdj[x, ] > 0))
                                                            cbind(Selidx, GrpAdj[x, Selidx])})
  #
  modified <- TRUE
  nb_pass_done <- 0

  while(modified & nb_pass_done != -1){
    cur_mod <- new_mod
    modified <- FALSE
    nb_pass_done <- nb_pass_done + 1
    rand_nodes <- sample(seq(length(node2com)), length(node2com), replace = FALSE)
    # print("***********************")

    for(u in rand_nodes){
      best_com <- node2com[u]
      best_increase <- 0
      #
      # print(node2com)
      weights2com <- data.frame(node= as.character(node2com[nbrs[[u]][, 1]]), sts= nbrs[[u]][, 2])
      weights2com <- aggregate(weights2com$sts, list(weights2com$node), FUN= sum)
      weights2com <- cbind(as.numeric(weights2com[, 1]), weights2com[, 2])
      #
      degrees <- degrees.v[u]
      Stot[best_com] <- Stot[best_com] - degrees
      #
      for(kk in sample(seq(dim(weights2com)[1]),  dim(weights2com)[1], replace = FALSE)){
        nbr_com <- weights2com[kk, 1]
        wt <- weights2com[kk, 2]
        #
        gain <- ((2*wt*norm) - (resolution*Stot[nbr_com] * degrees*2*norm^2))

        # print(c(kk, wt, gain))
        if(gain > best_increase){
          # print(c(nbr_com,gain))
          best_increase <- gain
          best_com <- nbr_com
        }
      }
      #
      Stot[best_com] <- Stot[best_com] + degrees
      #
      if(best_com != node2com[u]){
        #
        inner_partition[[node2com[u]]] <- as.vector(inner_partition[[node2com[u]]])[-match(as.vector(u), as.vector(inner_partition[[node2com[u]]]))]
        inner_partition[[best_com]] <- union(as.vector(inner_partition[[best_com]]), as.vector(u))
        modified <- TRUE
        node2com[u] <- best_com
      }
    }

    Grph_status$commu <- lapply(seq(length(Grph_status$commu)), function(zz) {
      if(length(inner_partition[[zz]] > 0)){
        unlist(lapply(seq(length(inner_partition[[zz]])), function(kk) sort(cur_Grph_status$commu[[inner_partition[[zz]][kk]]])))
      }else{
        integer(0)
      }
    }
    )

    new_mod <- LNModularity(Grph_status, resolution, norm)
    print(new_mod)
    #
    if((new_mod - cur_mod) <= threshold){
      inner_partition <- inner_partition[which(lengths(inner_partition)!=0)]
      clusterInfo <- Grph_status$commu[which(lengths(Grph_status$commu)!=0)]
      return(list(clusterInfo= clusterInfo, inner_partition= inner_partition, modified= modified, Modu = new_mod, GrpSta = Grph_status))
    }
  }
}

#'
#####****$$$$$$$$$$$$***Graph Partition**Second***$$$$$$$$$$$***#####
Lone_levelSecond <- function(G, resolution, threshold, Grph_status){
  #
  norm <- Grph_status$Norm
  #
  GrpAdj <- as.matrix(get.adjacency(G, attr = "weight", sparse = T))
  diag(GrpAdj) <- diag(GrpAdj)*2
  out_degree <- rowSums(GrpAdj)
  degrees.v <- out_degree
  print(sum(degrees.v))
  Stot <- degrees.v
  #
  vlist <- V(G)
  node2com <- as.vector(vlist)
  names(node2com) <- seq(length(node2com))
  inner_partition <- as.list(vlist)
  #
  cur_mod <- LNModularity(Grph_status, resolution, norm)
  cur_Grph_status <- Grph_status
  new_mod <- cur_mod
  #
  diag(GrpAdj) <- 0
  # nbrs <- lapply(1:length(node2com), function(x) cbind(as.numeric(which(GrpAdj[x, ] > 0)), GrpAdj[x, GrpAdj[x, ] > 0]))
  nbrs <- lapply(1:length(node2com), function(x) {Selidx <- which(GrpAdj[x, ] > 0)
                                                  cbind(Selidx, GrpAdj[x, Selidx])})
  #
  modified <- TRUE
  nb_pass_done <- 0

  while(modified & nb_pass_done != -1){
    cur_mod <- new_mod
    modified <- FALSE
    nb_pass_done <- nb_pass_done + 1
    rand_nodes <- sample(seq(length(node2com)), length(node2com), replace = FALSE)
    # print("***********************")

    for(u in rand_nodes){
      best_com <- node2com[u]
      best_increase <- 0
      #
      # print(node2com)
      weights2com <- data.frame(node= as.character(node2com[nbrs[[u]][, 1]]), sts= nbrs[[u]][, 2])
      weights2com <- aggregate(weights2com$sts, list(weights2com$node), FUN= sum)
      weights2com <- cbind(as.numeric(weights2com[, 1]), weights2com[, 2])
      #
      degrees <- degrees.v[u]
      Stot[best_com] <- Stot[best_com] - degrees
      #
      for(kk in sample(seq(dim(weights2com)[1]),  dim(weights2com)[1], replace = FALSE)){
        nbr_com <- weights2com[kk, 1]
        wt <- weights2com[kk, 2]
        #
        gain <- ((2*wt*norm) - (resolution*Stot[nbr_com] * degrees*2*norm^2))
        # print(c(kk, wt, gain))
        if(gain > best_increase){
          # print(c(nbr_com,gain))
          best_increase <- gain
          best_com <- nbr_com
        }
      }
      #
      Stot[best_com] <- Stot[best_com] + degrees
      #
      if(best_com != node2com[u]){
        #
        inner_partition[[node2com[u]]] <- as.vector(inner_partition[[node2com[u]]])[-match(as.vector(u), as.vector(inner_partition[[node2com[u]]]))]
        inner_partition[[best_com]] <- union(as.vector(inner_partition[[best_com]]), as.vector(u))
        modified <- TRUE
        node2com[u] <- best_com
      }
    }

    Grph_status$commu <- lapply(seq(length(Grph_status$commu)), function(zz) {
      if(length(inner_partition[[zz]] > 0)){
        unlist(lapply(seq(length(inner_partition[[zz]])), function(kk) sort(cur_Grph_status$commu[[inner_partition[[zz]][kk]]])))
      }else{
        integer(0)
      }
    }
    )

    new_mod <- LNModularity(Grph_status, resolution, norm)
    print(new_mod)
    #
    if((new_mod - cur_mod) <= threshold){
      # partition <- partition[which(lengths(partition)!=0)]
      inner_partition <- inner_partition[which(lengths(inner_partition)!=0)]
      clusterInfo <- Grph_status$commu[which(lengths(Grph_status$commu)!=0)]
      # return(list(partition=partition, inner_partition= inner_partition, modified= modified, Modu = new_mod))
      return(list(clusterInfo= clusterInfo, inner_partition= inner_partition, modified= modified, Modu = new_mod))
    }
  }
}

#'
######************Modularity Calculation***********************######
LNModularity <- function(Grph_status, resolution, norm){
  #
  G <- Grph_status$Grp
  commu <- Grph_status$commu
  out_degree <- Grph_status$ODeg
  in_degree <- Grph_status$IDeg
  norm <- Grph_status$Norm
  #
  if(!is.list(commu)){
    stop("Graph partition is not a list")
  }
  #
  if(vcount(G) != length(unique(unlist(commu)))){
    stop("Partition is not the Graph partition")
  }
  #
  modulSum <- 0
  for(ii in seq(length(commu))){
    community <- as.vector(commu[[ii]])
    IndG <- induced_subgraph(G, community)
    #
    out_degree_sum <- sum(out_degree[community])
    in_degree_sum <- sum(in_degree[community])
    #
    modulSum <- modulSum + ((sum(strength(IndG))*norm) - (resolution * out_degree_sum * in_degree_sum * (norm^2)))
  }

  return(modulSum)
}

#'
##############Generate Updated Graph##############
Lgen_graph <- function(G, partition){
  # """Generate a new graph based on the partitions of a given graph"""
  Cno <- length(partition)
  EdgInfo <- as_edgelist(G)
  EdgInfo <- cbind(EdgInfo, E(G)$weight, matrix(data = NA, nrow = dim(EdgInfo)[1], ncol = 2))
  Cinfo <- cbind(unlist(partition), unlist(lapply(1:Cno, function(i) rep(i, length(partition[[i]])))))
  for(kk in seq(dim(Cinfo)[1])){
    EdgInfo[which(Cinfo[kk, 1] == EdgInfo[, 1]), 4] <- Cinfo[kk, 2]
    EdgInfo[which(Cinfo[kk, 1] == EdgInfo[, 2]), 5] <- Cinfo[kk, 2]
  }
  df <- data.frame(EdgInfo[, c(3:5)]) %>% group_by(X2, X3) %>% summarise(EW= sum(as.numeric(X1)))

  H <- graph_from_edgelist(cbind(as.numeric(df$X2), as.numeric(df$X3)), directed = FALSE)
  E(H)$weight <- df$EW
  return(H)
}
