#' @title
#' ELVAR algorithm detects communities from a node-attributed graph network
#'
#' @aliases Eva_partitions
#'
#' @description
#' The code is a generalized Louvain clustering method that we call ELVAR (Extended Louvain Algorithm),
#' which is an R-implementation of the EVA algorithm. In ELVAR algorithm, node attributes information
#' along with graph topology enhance to find more relevant clusters specially when
#' clusters are not well-separated in high dimensional space.
#'
## Use Eva_partitions() function to run the code...
#'
## Input parameters......
#' @param G
#' A node-attributed igraph object
#' @param resolution
#' Tune the size of the communities. Default is 1. Community size decreases with the increment of resolution value
#' @param threshold
#' If change in the quality function (score) is less than the threshold, iteration loop will be terminated
#' @param alpha
#' Tuning parameter that makes balance between purity and modularity during optimization
#' @param order
#' Defines numerical order among the characters present within an attribute
#' @param Vattr.name
#' User-select node attributes that will be considered for Eva algorithm
#'
## Output parameters......
#' @return
#' It is a list consists of
#' CommunityMembers:: is a vector where names are the nodes and the values are the communities belong to the nodes,
#' ModularityVal:: Modularity value, and
#' PurityVal:: Purity value
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
Eva_partitions <- function(G, resolution= 1, threshold= 0.0000001, alpha=0.5, order=NULL, Vattr.name){

  ####

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
  }
  # Delete vertex attributes that are not defined by the user
  remvAttri <- setdiff(vertex_attr_names(G), Vattr.name)
  if(length(remvAttri)!=0){
    for(rmv in remvAttri){
      G <- delete_vertex_attr(G, rmv)
    }
  }
  # Run Eva optimization
  OneLevelSumm <- one_levelFirst(G, resolution, threshold, alpha, order)
  #
  mod <- OneLevelSumm$Modu
  pur <- OneLevelSumm$Purit
  CommInfo <- OneLevelSumm$clusterInfo
  Grph_status <- OneLevelSumm$GrpSta
  Grph_status$commu <- CommInfo
  #
  current_grp <- gen_graph(G, CommInfo)
  #
  while(TRUE){
    OneLevelSumm <- one_levelSecond(current_grp, resolution, threshold, Grph_status, alpha, order)
    #
    new_mod <- OneLevelSumm$Modu
    new_pur <- OneLevelSumm$Purit
    CommInfo <- OneLevelSumm$clusterInfo
    Grph_status$commu <- CommInfo

    score = alpha * (new_pur - pur) + (1 - alpha) * (new_mod - mod)
    print(score)

    if(score <= threshold){

      members <- rep(1:length(CommInfo), lengths(CommInfo))
      names(members) <- unlist(CommInfo)
      members <- members[sort(as.numeric(names(members)), index.return= TRUE)$ix]
      if(!is.null(vertexN.v)){
        names(members) <- vertexN.v
        }
      return(list(CommunityMembers= members, ModularityVal= new_mod, PurityVal= new_pur))
    }
    #
    mod <- new_mod
    pur <- new_pur
    current_grp <- gen_graph(G, CommInfo)
   }
}

#### Auxiliary functions
#'
#####****$$$$$$$$$$$$***Graph Partition***First**$$$$$$$$$$$***#####
one_levelFirst <- function(G, resolution, threshold, alpha, order){
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
  cur_mod <- NModularity(Grph_status, resolution, norm)
  new_mod <- cur_mod
  #
  ClusatriInfo <- VrtAttriSumm(Grph_status)
  cur_purity <- overall_purity(ClusatriInfo)
  new_purity <- cur_purity

  # print(c(new_mod, new_purity))
  # print("%%%%%%%%%%%%")
  cur_Grph_status <- Grph_status
  #
  diag(GrpAdj) <- 0
  nbrs <- lapply(1:length(node2com), function(x) {Selidx <- as.numeric(which(GrpAdj[x, ] > 0))
                                                  cbind(Selidx, GrpAdj[x, Selidx])})
  #
  modified <- TRUE
  nb_pass_done <- 0

  while(modified & nb_pass_done != -1){
    cur_mod <- new_mod
    cur_purity <- new_purity

    modified <- FALSE
    nb_pass_done <- nb_pass_done + 1
    rand_nodes <- sample(seq(length(node2com)), length(node2com), replace = FALSE)
    # print("***********************")

    for(u in rand_nodes){
      best_mod <- 0
      best_com <- node2com[u]
      best_size_incr <- 0
      #
      # print(node2com)
      weights2com <- data.frame(node= as.character(node2com[nbrs[[u]][, 1]]), sts= nbrs[[u]][, 2])
      weights2com <- aggregate(weights2com$sts, list(weights2com$node), FUN= sum)
      weights2com <- cbind(as.numeric(weights2com[, 1]), weights2com[, 2])
      degrees <- degrees.v[u]
      Stot[best_com] <- Stot[best_com] - degrees

      #
      for(kk in sample(seq(dim(weights2com)[1]),  dim(weights2com)[1], replace = FALSE)){
        # gain by node attributes
        initital_labels <- ClusatriInfo[[best_com]]

        nbr_com <- weights2com[kk, 1]
        wt <- weights2com[kk, 2]
        #
        incr_labels <- ClusatriInfo[[nbr_com]]
        delPurity <- delta_purity_size(initital_labels, incr_labels, order)
        incr_attr <- delPurity[1]
        incr_size <- delPurity[2]
        #
        gain <- ((2*wt*norm) - (resolution*Stot[nbr_com] * degrees*2*norm^2))*10^3
        #
        total_incr <- alpha*incr_attr + (1 - alpha)*gain    # Quality function
        # print(c(nbr_com, wt, incr_attr, gain, total_incr))

        if(total_incr > best_mod | (total_incr == best_mod & incr_size > best_size_incr)){
          best_mod <- total_incr
          best_size_incr <- incr_size
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

    new_mod <- NModularity(Grph_status, resolution, norm)
    #
    ClusatriInfo <- VrtAttriSumm(Grph_status)
    new_purity <- overall_purity(ClusatriInfo)
    #
    score <- alpha * (new_purity - cur_purity) + (1 - alpha) * (new_mod - cur_mod)
    print(c(score, new_mod, new_purity))

    #
    if(score <= threshold){
      # partition <- partition[which(lengths(partition)!=0)]
      inner_partition <- inner_partition[which(lengths(inner_partition)!=0)]
      clusterInfo <- Grph_status$commu[which(lengths(Grph_status$commu)!=0)]
      # return(list(partition=partition, inner_partition= inner_partition, modified= modified, Modu = new_mod))
      return(list(clusterInfo= clusterInfo, inner_partition= inner_partition, modified= modified, Modu = new_mod,
                  Purit= new_purity, GrpSta = Grph_status))
    }
  }
}

#'
#####****$$$$$$$$$$$$***Graph Partition***Second**$$$$$$$$$$$***#####
one_levelSecond <- function(G, resolution, threshold, Grph_status, alpha, order){
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
  cur_mod <- NModularity(Grph_status, resolution, norm)
  new_mod <- cur_mod
  #
  ClusatriInfo <- VrtAttriSumm(Grph_status)
  cur_purity <- overall_purity(ClusatriInfo)
  new_purity <- cur_purity

  # print(c(new_mod, new_purity))
  # print("%%%%%%%%%%%%")
  cur_Grph_status <- Grph_status
  #
  diag(GrpAdj) <- 0
  nbrs <- lapply(1:length(node2com), function(x) {Selidx <- which(GrpAdj[x, ] > 0)
                                                  cbind(Selidx, GrpAdj[x, Selidx])})
  #
  modified <- TRUE
  nb_pass_done <- 0

  while(modified & nb_pass_done != -1){
    cur_mod <- new_mod
    cur_purity <- new_purity

    modified <- FALSE
    nb_pass_done <- nb_pass_done + 1
    rand_nodes <- sample(seq(length(node2com)), length(node2com), replace = FALSE)
    # print("***********************")

    for(u in rand_nodes){
      best_mod <- 0
      best_com <- node2com[u]
      best_size_incr <- 0
      #
      # print(node2com)
      weights2com <- data.frame(node= as.character(node2com[nbrs[[u]][, 1]]), sts= nbrs[[u]][, 2])
      weights2com <- aggregate(weights2com$sts, list(weights2com$node), FUN= sum)
      weights2com <- cbind(as.numeric(weights2com[, 1]), weights2com[, 2])
      degrees <- degrees.v[u]
      Stot[best_com] <- Stot[best_com] - degrees

      #
      for(kk in sample(seq(dim(weights2com)[1]),  dim(weights2com)[1], replace = FALSE)){
        # gain by node attributes
        initital_labels <- ClusatriInfo[[best_com]]

        nbr_com <- weights2com[kk, 1]
        wt <- weights2com[kk, 2]
        #
        incr_labels <- ClusatriInfo[[nbr_com]]
        delPurity <- delta_purity_size(initital_labels, incr_labels, order)
        incr_attr <- delPurity[1]
        incr_size <- delPurity[2]
        #
        gain <- ((2*wt*norm) - (resolution*Stot[nbr_com] * degrees*2*norm^2))*10^3
        #
        total_incr <- alpha*incr_attr + (1 - alpha)*gain    # Quality function
        # print(c(kk, wt, incr_attr, gain))
        if(total_incr > best_mod | (total_incr == best_mod & incr_size > best_size_incr)){
          best_mod <- total_incr
          best_size_incr <- incr_size
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

    new_mod <- NModularity(Grph_status, resolution, norm)
    #
    ClusatriInfo <- VrtAttriSumm(Grph_status)
    new_purity <- overall_purity(ClusatriInfo)
    #
    score <- alpha * (new_purity - cur_purity) + (1 - alpha) * (new_mod - cur_mod)
    print(c(score, new_mod, new_purity))

    #
    if(score <= threshold){
      # partition <- partition[which(lengths(partition)!=0)]
      inner_partition <- inner_partition[which(lengths(inner_partition)!=0)]
      clusterInfo <- Grph_status$commu[which(lengths(Grph_status$commu)!=0)]
      # return(list(partition=partition, inner_partition= inner_partition, modified= modified, Modu = new_mod))
      return(list(clusterInfo= clusterInfo, inner_partition= inner_partition, modified= modified, Modu = new_mod, Purit= new_purity))
    }
  }
}

#'
######************Modularity Calculation***********************######
NModularity <- function(Grph_status, resolution, norm){
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
gen_graph <- function(G, partition){
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

#'
#######################Local purity calculation#################
delta_purity_size <- function(original_atr1, new_atr2, order){

  total_original <-  sum(unlist(original_atr1))/length(original_atr1)
  total_nodes <- total_original

  dumping <- 1
  if(!is.null(order)){
    for(atr in names(original_atr1)){
      most_freq_label <- names(original_atr1[[atr]])[match(max(original_atr1[[atr]]), original_atr1[[atr]])]
      most_freq_label_new <- names(new_atr2[[atr]])[match(max(new_atr2[[atr]]), new_atr2[[atr]])]
      hierarchy <- order[[atr]]
      d <- 1 - as.numeric(abs(hierarchy[most_freq_label] - hierarchy[most_freq_label_new]))/length(hierarchy)
      dumping <- dumping*d
    }
  }

  com_pur <- c()
  for(jj in seq(length(original_atr1))){
    atri.v <- original_atr1[[jj]]
    score <- atri.v/sum(atri.v)
    com_pur <- append(com_pur, max(score))
  }
  purity_original <- prod(com_pur, na.rm = TRUE)

  # computing original purity
  total_nodes <- total_nodes + sum(unlist(new_atr2))/length(new_atr2)

  updated <- original_atr1
  for(ata in names(new_atr2)){
    if(ata %in% names(updated)){
      for(tv in names(new_atr2[[ata]])){
        if(tv %in% names(updated[[ata]])){
          updated[[ata]][tv] <- updated[[ata]][tv] + new_atr2[[ata]][tv]
        }else{
          updated[[ata]][tv] <- new_atr2[[ata]][tv]
        }
      }
    }else{
      updated[[ata]] <- new_atr2[[ata]]
    }
  }

  com_pur <- c()
  for(jj in seq(length(updated))){
    atri.v <- updated[[jj]]
    score <- atri.v/sum(atri.v)
    com_pur <- append(com_pur, max(score))
  }
  purity_overall <- prod(com_pur, na.rm = TRUE)*dumping
  #
  increment <- purity_overall - purity_original
  delta_size <- (total_nodes - total_original)/total_original

  return(c(increment, delta_size))
}

#'
##############Global purity calculation#############
overall_purity <- function(ClusatriInfo){
  purities <- c()
  for(kk in seq(length(ClusatriInfo))){
    com_pur <- c()
    if(lengths(ClusatriInfo[[kk]])[1]==0){
      com_pur <- rep(NA, length(ClusatriInfo[[kk]]))
    }else{
      for(jj in seq(length(ClusatriInfo[[kk]]))){
        atri.v <- ClusatriInfo[[kk]][[jj]]
        score <- atri.v/sum(atri.v)
        com_pur <- append(com_pur, max(score))
      }
    }
    purities <- append(purities, prod(com_pur, na.rm = TRUE))
  }
  return(mean(purities))
}

#'
############Update community attribute summary##########
VrtAttriSumm <- function(Grph_status){
  G <- Grph_status$Grp
  partition <- Grph_status$commu
  ClusatriInfo <- list()
  for(jj in seq(length(partition))){
    Eachclustatri <- vertex.attributes(G, partition[[jj]])
    atr <- names(Eachclustatri)
    ClusatriInfo[[jj]] <- lapply(atr, function(atr) summary(factor(Eachclustatri[[atr]]), maxsum = "all"))
    names(ClusatriInfo[[jj]]) <- atr
  }
  return(ClusatriInfo)
}
