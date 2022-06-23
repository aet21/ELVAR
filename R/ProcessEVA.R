#' @title 
#' Processes the output of the EVA clustering algorithm
#' 
#' @aliases ProcessEVA
#'  
#' @description 
#' This function takes as input the result of an EVA-clustering plus the metacell information from the seurat
#' scRNA-Seq data object, to identify for each attribute value the communities enriched for cells with that
#' attribute value. This is done using a binomial test. Once the enriched clusters for each attribute value
#' are found, cells from these clusters and with the given attribute value are merged.
#' 
#' @param eva.o
#' The output from running the EVA-clustering algorithm.
#' 
#' @param seu.o
#' The Seurat scRNA-Seq data object containing the metacell information.
#'
#' @param attrName
#' The name of the attribute to use for testing enrichment of clusters. Must be one of the names
#' in the metadata slot of the Seurat object.
#'
#' @param sigth
#' The significance threshold on the binomial-test P-values for declaring enrichment.
#' If not specified (default), the function will use a Bonferroni correction.
#' 
#' @return A list with four elements
#' 
#' @return pvalEnr
#' The enrichment P-value matrix defined over clusters (rows) and attribute values (columns).
#'
#' @return sigth
#' The significance threshold used to declare enriched clusters.
#' 
#' @return sigClust
#' A list of enriched cluster identifiers, each list entry representing one attribute value.
#' 
#' @return cellsMrg
#' A list of cell indices indicating which cells of a given attribute value are in enriched clusters
#' for that attribute value. Each entry in the list represents an attribute value.
#' 
#' 
#' @references 
#' Qi Luo, Alok K Maity, Andrew E Teschendorff
#' \emph{Distance Covariance Entropy reveals primed states and bifurcation dynamics in single-cell RNA-Seq data}
#' Submitted.
#' 
#' 
#' @examples 
#'
#' 
#' @export
#'

ProcessEVA <- function(eva.o,seu.o,attrName,sigth=NULL){

   tmp.idx <- which(names(seu.o@meta.data)==attrName);
   if(length(tmp.idx)==0){
        print("Attribute name not found! Please recheck!");
        break;
   }
   else if (length(tmp.idx)>0){
       tabClust.m <- table(eva.o$CommunityMember,seu.o@meta.data[[tmp.idx]]);
   }
   
  nCellAttrVal.v <- apply(tabClust.m,2,sum)
  nCellClust.v <- apply(tabClust.m,1,sum)
  nCells <- sum(nCellAttrVal.v);

  probAttr.v <- nCellAttrVal.v/nCells;
  pvalEnrAttrVal.m <- tabClust.m;
  for(attr in 1:length(nCellAttrVal.v)){
    for(cl in 1:nrow(tabClust.m)){
      if(tabClust.m[cl,attr] < nCellClust.v[cl]){
          pvalEnrAttrVal.m[cl,attr] <- pbinom(tabClust.m[cl,attr],size=nCellClust.v[cl],prob=probAttr.v[attr],lower.tail=FALSE);
      }
      else if (tabClust.m[cl,attr]==nCellClust.v[cl]){
          pvalEnrAttrVal.m[cl,attr] <- probAttr.v[attr]^nCellClust.v[cl];
      }
    }
  }

  if(is.null(sigth)){
   sigth <- 0.05/(length(nCellClust.v)*length(nCellAttrVal.v));
  }
    
  sigClustAttrVal.lv <- list();
  for(attr in 1:length(nCellAttrVal.v)){
   sig.idx <- which(pvalEnrAttrVal.m[,attr] < sigth);
   sigClustAttrVal.lv[[attr]] <- names(sig.idx);
  }

 nClustAttrVal.v <- unlist(lapply(sigClustAttrVal.lv,length));
    
 ### Now merge cells from different enriched clusters per attribute value
 attrVal.v <- colnames(tabClust.m);
 cellsMrg.li <- list();
 for(attr in 1:length(nCellAttrVal.v)){
   cellsMrg.li[[attr]] <- intersect(which(eva.o$CommunityMembers %in% sigClustAttrVal.lv[[attr]]),which(seu.o@meta.data[[tmp.idx]]==attrVal.v[attr]));
 }

  return(list(pvalEnr=pvalEnrAttrVal.m,sigth=sigth,sigClust=sigClustAttrVal.lv,cellsMrg=cellsMrg.li));
}

