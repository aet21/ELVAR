#' @title 
#' Performs DA-analysis
#' 
#' @aliases DoDA
#'  
#' @description 
#' This function performs differential abundance (DA) testing of a cell-state variable specified
#' by \code{varDA} across the cell-attribute used in the previous EVA-clustering step. Replicate
#' information is specified via the \code{varREP} variable. The cell annotations must be present in
#' the metadata slots of the Seurat scRNA-Seq data object \code{seu.o}. The other inputs to this function
#' include the output of one EVA-clustering run, as well as the output from \code{ProcessEVA}. The function
#' then returns the negative binomial regression statistics to determine if the cell-state fractions vary
#' significantly as a function of attribute value.
#' 
#' @param eva.o
#' The output from running the EVA-clustering algorithm.
#' 
#' @param seu.o
#' The Seurat scRNA-Seq data object containing the metacell information.
#'
#' @param sigcl.o
#' The output from running the \code{ProcessEVA} function.
#'
#' @param varDA
#' A character specifying the name of the cell-state for which DA is being
#' tested for. 
#'
#' @param varREP
#' A character specifying the name of replicate ID variable.
#' 
#' @return A list with two elements
#' 
#' @return stat
#' A matrix of statistics from the negative binomial regression, with each row
#' representing a cell-state and columns labeling the z-statistic and P-value.
#'
#' @return df
#' The data frame used for the negative binomial regression. Each row labels one
#' attribute-replicate pair. Columns label the attribute value, the two states for
#' which DA is to be tested across the attribute, and the normalization factor.
#'
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
#' @import MASS
#' 
#' 
#' @export
#'

DoDA <- function(eva.o,seu.o,sigcl.o,varDA,varREP){

   varDA.idx <- which(names(seu.o@meta.data)==varDA);
   if(length(varDA.idx)==0){
        print("Variable for DA not found! Please recheck!");
        break;
   }
  
   varREP.idx <- which(names(seu.o@meta.data)==varREP);
   if(length(varREP.idx)==0){
        print("Variable defining sample replicates not found! Please recheck!");
        break;
   }

   nClustAttrVal.v <- unlist(lapply(sigcl.o$sigClust,length));
   selAttrVal.idx <- which(nClustAttrVal.v>0);
    
   npm <- 0;
   for(attr in selAttrVal.idx){
     npm <- npm + length(unique(seu.o@meta.data[[varREP.idx]][sigcl.o$cellsMrg[[attr]]]));
   }

 
   states.v <- unique(seu.o@meta.data[[varDA.idx]]);
   stateN.v <- paste("CS",states.v,sep="");
   ncol <- 2+length(states.v);
   tmpDFeva.m <- matrix(NA,nrow=npm,ncol=ncol);
   colnames(tmpDFeva.m) <- c("Attr",stateN.v,"NF");

   rowIDX <- 1;
   for(attr in selAttrVal.idx){
      st1.idx <- which(seu.o@meta.data[[varDA.idx]][sigcl.o$cellsMrg[[attr]]]==states.v[1]);
      st2.idx <- which(seu.o@meta.data[[varDA.idx]][sigcl.o$cellsMrg[[attr]]]==states.v[2]);      
      unqRepID.v <- unique(seu.o@meta.data[[varREP.idx]][sigcl.o$cellsMrg[[attr]]]);
      for(r in 1:length(unqRepID.v)){
       tmp.idx <- which(seu.o@meta.data[[varREP.idx]][sigcl.o$cellsMrg[[attr]]]==unqRepID.v[r]);
       tmpDFeva.m[rowIDX,1] <- attr;
       tmpDFeva.m[rowIDX,2] <- length(intersect(tmp.idx,st1.idx));
       tmpDFeva.m[rowIDX,3] <- length(intersect(tmp.idx,st2.idx));
       tmpDFeva.m[rowIDX,4] <- log(length(tmp.idx));
       rowIDX <- rowIDX +1;
      }
   }

   statNBR.m <- matrix(nrow=length(stateN.v),ncol=2);
   rownames(statNBR.m) <- stateN.v;
   colnames(statNBR.m) <- c("z","P");    
   for(i in 1:length(states.v)){
    nbr.o <- glm.nb( tmpDFeva.m[,i+1] ~ Attr + NF, data=as.data.frame(tmpDFeva.m)) ; ### default link is log (what we want)
    statNBR.m[i,1:2] <- summary(nbr.o)$coeff[2,3:4];
   }
 
   return(list(stat=statNBR.m,df=tmpDFeva.m));
}

