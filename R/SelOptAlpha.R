#' @title
#' Selects optimal purity parameter value
#' 
#' @aliases SelOptAlpha
#'
#' @description
#' The aim of this function is to facilitate choosing an optimal purity parameter (alpha) value.
#' As input it requires the number of inferred clusters, the purity and modularity values for a range
#' of different alpha values and for a reasonably large number of runs (typically we recommend
#' 100 runs), as well as the number of clusters inferred with the ordinary Louvain algorithm.
#' If a non-sequential stochastic version of Louvain is used, then this could be an average
#' value over multiple runs.
#'
#' @param nC.m
#' A matrix of inferred cluster numbers, with rows labeling runs and columns labeling
#' the alpha-values.
#' 
#' @param nc0
#' A scalar giving the number of inferred clusters using the Louvain algorithm.
#'
#' @param pur.m
#' A matrix of estimated purity values, with rows labeling the runs and columns labeling
#' the alpha-values.
#'
#' @param mod.m
#' A matrix of estimated modularity values, with rows labeling the runs and columns labeling
#' the alpha-values.
#'
#' @param alpha.v
#' A vector of purity parameter values. Should be an increasing sequence between 0 and 1,
#' excluding these values as these are not interesting. By default, this is a sequence from
#' 0.1 to 0.9. We don't advise changing the extreme values, but user can change the granularity
#' or step-size.
#'
#' @param thP
#' A scalar between 0 and 1, which quantifies by how much purity at the optimal alpha-value
#' is allowed to deviate from the maximum attained over the whole range of alpha-values
#' considered. By default this value is 0.75. A value of 1 would mean no deviation from the maximum.
#'
#' @param thM
#' A scalar between 0 and 1, which quantifies by how much the modularity at the optimal alpha-value
#' is allowed to deviate from the maximum attained over the whole range of alpha-values considered.
#' By default, this values is 0.75. A value of 1 would mean no deviation from the maximum.
#'
#' @param qC
#' A scalar between 0 and 1 labeling the quantile of the EVA cluster number distribution to consider.
#' By default, this value is 0.95, which means that the upper 95 percent quantile of the EVA cluster
#' number distribution is used when comparing to \code{nc0}.
#'
#' @param thC
#' A scalar quantifying by how much higher the EVA cluster number obtained via \code{qC} should
#' be in relation to \code{nc0}. By default, this value is 1.5, which means that we require the
#' upper 95 percent quantile of the EVA cluster number distribution to be at least 1.5 higher than
#' \code{nc0}.
#'
#' @return A list with three elements
#'
#' @return aM
#' The largest allowed alpha-value satisfying the modularity constraint.
#'
#' @return aP
#' The smallest allowed alpha-value satisfying the purity constraint.
#'
#' @return aC
#' The smallest allowed alpha-value satisfying the cluster number constraint.
#'
#'
#' @references
#' Alok K Maity, Andrew E Teschendorff
#' \emph{Cell-Attribute aware community detection improves differential abundance testing from single-cell RNA-Seq data}. Submitted.
#'
#' @examples
#'
#' @export


SelOptAlpha <- function(nC.m,nc0,pur.m,mod.m,alpha.v=seq(0.1,0.9,0.1),thP=0.75,thM=0.75,qC=0.95,thC=1.5){

### for purity select smallest alpha for which purity is within 75% of their maximum value
### for modularity select largest alpha for which modularity is within 75% of their maximum value    
    avP.v <- colMeans(pur.m);
    avM.v <- colMeans(mod.m);
    maxP <- max(avP.v);
    maxM <- max(avM.v);
    fP.v <- avP.v/maxP;
    fM.v <- avM.v/maxM;

    aP <- alpha.v[min(which(fP.v > thP))];
    aM <- alpha.v[max(which(fM.v > thM))];    

    rC.m <- nC.m/nc0;
    qC.v <- apply(rC.m,2,quantile,qC);
    aC <- alpha.v[min(which(qC.v >= thC))];

    
    return(list(aM=aM,aP=aP,aC=aC));
}
