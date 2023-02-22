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
#' then returns the negative binomial regression statistics (Wald and LRT) to determine if the cell-state fractions vary
#' significantly as a function of attribute value.
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
#' A matrix of statistics (Wald and LRT) from the negative binomial regression, with each row
#' representing a cell-state and columns labeling the z-statistic and P-value corresponding to
#' the statistical tests.
#'
#' @return df
#' The data frame used for the negative binomial regression. Each row labels one
#' attribute-replicate pair. Columns label the attribute value, the two states for
#' which DA is to be tested across the attribute, and the normalization factor.
#'
#'
#'
#' @references
#' Alok K Maity, Andrew E Teschendorff
#' \emph{Cell-Attribute aware community detection improves differential abundance testing from single-cell RNA-Seq data}
#' Submitted.
#'
#'
#' @examples
#'
#' @import MASS
#'
#'
#' @export


DoDA <- function(seu.o,sigcl.o,varDA,varREP){

  varDA.idx <- which(names(seu.o@meta.data) == varDA)
  if (length(varDA.idx) == 0) {
    stop("Variable for DA not found! Please recheck!")
  }
  varREP.idx <- which(names(seu.o@meta.data) == varREP)
  if (length(varREP.idx) == 0) {
    stop("Variable defining sample replicates not found! Please recheck!")
  }
  nClustAttrVal.v <- unlist(lapply(sigcl.o$sigClust, length))
  selAttrVal.idx <- which(nClustAttrVal.v > 0)
  npm <- 0
  for (attr in selAttrVal.idx) {
    npm <- npm + length(unique(seu.o@meta.data[[varREP.idx]][sigcl.o$cellsMrg[[attr]]]))
  }
  states.v <- unique(seu.o@meta.data[[varDA.idx]])
  stateN.v <- paste("CS", states.v, sep = "")
  ncol <- 2 + length(states.v)
  tmpDFeva.m <- matrix(NA, nrow = npm, ncol = ncol)
  colnames(tmpDFeva.m) <- c("Attr", stateN.v, "NF")
  rowIDX <- 1
  for (attr in selAttrVal.idx) {
    st.idx <- lapply(states.v, function(states.v) which(seu.o@meta.data[[varDA.idx]][sigcl.o$cellsMrg[[attr]]]==states.v))

    unqRepID.v <- unique(seu.o@meta.data[[varREP.idx]][sigcl.o$cellsMrg[[attr]]])

    for (r in 1:length(unqRepID.v)) {
      tmp.idx <- which(seu.o@meta.data[[varREP.idx]][sigcl.o$cellsMrg[[attr]]] ==
                         unqRepID.v[r])
      tmpDFeva.m[rowIDX, 1] <- attr
      tmpDFeva.m[rowIDX,(c(1:length(st.idx))+1)] <- sapply(st.idx, function(st.idx) length(intersect(tmp.idx,st.idx)))

      tmpDFeva.m[rowIDX, (length(st.idx)+2)] <- log(length(tmp.idx));
      rowIDX <- rowIDX + 1
    }
  }
  #
  if(length(summary(factor(tmpDFeva.m[, 1]))) < 2){
    stop("Enriched Communities come from only one disease state or age group!");
  }
  #

  # Performed DA-analysis by Negative Binomial Regression and used Wald & LRT statistics
  statNBR.m <- matrix(nrow = length(stateN.v), ncol = 4)
  rownames(statNBR.m) <- stateN.v
  colnames(statNBR.m) <- c("Wald(z)", "Wald(P)", "LRT(z)", "LRT(P)")

  for(i in c(1:length(states.v))){
    nbr.o <- glm.nb(tmpDFeva.m[, i + 1] ~ Attr + NF, data = as.data.frame(tmpDFeva.m)) ### default link is log (what we want)

    chisq.o <- pchisq(abs(nbr.o$twologlik),df=1,lower.tail=FALSE)
    z_val <- qnorm(chisq.o, 0, 1, lower.tail = (sign(nbr.o$coef[2])<0))

    statNBR.m[i, ] <- c(summary(nbr.o)$coeff[2, 3:4], z_val, chisq.o)
  }

   return(list(stat=statNBR.m,df=tmpDFeva.m));
}
