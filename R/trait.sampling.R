#'Simulated sampling for trait value measurement
#'
#'Randomly selects individuals for trait value measurement and gives back raw measured
#'traits or their means
#'
#'It simulates the real world situation that not all individuals are collected for trait measurement.
#'If \code{ITV==FALSE}, all individuals belonging to the species are pooled, and then \code{n}
#' randomly selected individuals are measured. If \code{ITV==TRUE}, \code{n} individuals are
#' measured in each (sub)community, where the species occur. If the occurring individuals are
#' less than \code{n}, all individuals are measured.
#'
#'If \code{aggregate==TRUE}, meta-community or subcommunity level means are calculated, otherwise
#'raw measurements are returned.
#'
#'@param x    community and trait data matrix produced by \code{\link{comm.simul}}
#'            function
#'@param ITV  If \code{TRUE} each subcommunity are sampled separately, otherwise
#'            the meta-community level sampling was done
#'@param n    Number of sampled individuals
#'@param aggregate    If \code{TRUE} mean trait values are returned, otherwise the raw values
#'                    of sampled individuals
#'@return  data.frame with fields: \code{species}, \code{site} (only if \code{ITV=TRUE}),
#'        \code{trait.a}, \code{trait.b}, \code{trait.c} (raw values or means depending on parameter \code{aggregate})
#'
#'@export
#'
#'@examples
#' x<-comm.simul(S=20, J=30)
#' str(x)
#'
#'w<-trait.sampling(x$final.community)
#'w
#'
#'w<-trait.sampling(x$final.community,ITV=TRUE,aggregate=TRUE)
#'str(w)
#'
#'


trait.sampling<-function(x,ITV=FALSE,aggregate=TRUE,n=5)
{
  x2<-vector()
  sp_names<-levels(as.factor(x[,"species"]))
  if (ITV)
    {
    site_names<-levels(as.factor(x[,"site"]))
    for (i in sp_names)
      for (j in site_names)
        {
        w<-x[(x[,"species"]==i) & (x[,"site"]==j),]
        if (nrow(w)>n) w<-w[sample(1:nrow(w),n),]
        x2<-rbind(x2,w)
        }
    if (aggregate)  x2<-aggregate(x2[,3:5],by=list(Site=x2[,1],Species=x2[,2]),
                                  mean)
    }
  else
    {
    for (i in sp_names)
      {
      w<-x[x[,"species"]==i,]
      if (nrow(w)>n) w<-w[sample(1:nrow(w),n),]
      x2<-rbind(x2,w)
      }
    x2<-x2[,-1]
    if (aggregate)  x2<-aggregate(x2[,2:4],by=list(Species=x2[,1]),mean)
    }
rownames(x2)<-c()
return(x2)
}
