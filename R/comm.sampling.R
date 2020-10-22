#' Converting simulation results into site-by-species matrix
#'
#'Converts simulation result into site-by-species matrix of abundances,
#'and optionally in the same step simulates random sampling with fixed number of individuals.
#'
#'If \code{type=="full"} it simply converts simulation results from long to wide format.
#'If \code{type=="random"} it randomly selects \code{size} individuals in each (sub)community
#'and abundances in these samples are converted into site-by-species matrix format.
#'
#'@param x    community and trait data matrix produced by \code{\link{comm.simul}}
#'            function
#'@param type Type of sampling. If \code{type=="full"} sample size equals to community size.
#'            It simply converts \code{x} into a site-by-species matrix
#'            If \code{type=="random"}, it applies random sampling by calling (\code{\link[vegan]{rrarefy}}) function
#'@param size  Number of individuals in the random samples. It should be smaller
#'             than number of individuals in simulated (sub)communities. Otherwise,
#'             \code{x} is converted into a site-by-species matrix without (re)sampling
#'@return A site-by-species matrix containing abundances.
#'
#'@export
#'
#'@examples
#' x<-comm.simul(S=20, J=30)
#' str(x$final.community)
#'
#'w<-comm.sampling(x$final.community,type="full")
#'str(w)
#'
#'w.rarefied<-comm.sampling(x$final.community,type="random",size=10)
#'rowSums(w)
#'rowSums(w.rarefied)
#'


comm.sampling<-function(x,type=c("full","random"),size)
{
  type <- match.arg(type)
  x2<-as.matrix(table(x[,"site"],x[,"species"]))
  if (type=="random")
    {
    if (size>=min(rowSums(x2))) warning("Sample size should be smaller than actual community size! Full sampling will be applied.")
    else x2<-vegan::rrarefy(x2,size)
  }
 return(x2)
}
