\name{competition.kernel}
\alias{competition.kernel}
\title{Competition kernels}

\arguments{
\item{trait.values}{Values of trait related to resource use}

\item{...}{Additional parameters}
}
\value{
Square matrix of pairwise competition coefficients
}
\description{
User defined functions for calculating pairwise competition coefficients
}
\details{
User can defined any specific from of competition.
Pairwise competition between species/individuals should depend on their
trait values related to resource use. Vector of these trait values has to be
the first parameter of the function, and any further parameters are allowed.
The output has to be a square matrix of pairwise competition coefficients.


Competition kernels available in the package: \cr
\code{\link{asymmetric.competition.kernel}} \cr
\code{\link{Gaussian.competition.kernel}} \cr
}
