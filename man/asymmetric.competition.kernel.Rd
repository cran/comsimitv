% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asymmetric.competition.kernel.R
\name{asymmetric.competition.kernel}
\alias{asymmetric.competition.kernel}
\title{Asymmetric competition kernels}
\usage{
asymmetric.competition.kernel(
  trait.values,
  trait.compet = "trait.b",
  ac.type = c("Kisdi", "Nattrass"),
  sigma.b = 0.03,
  ac.C = 1,
  ac.v = 1,
  ...
)
}
\arguments{
\item{trait.values}{Dataframe of all traits}

\item{trait.compet}{Name of trait related to resource use}

\item{ac.type}{Type of the function (see \code{vignette("competition")})}

\item{sigma.b}{steepness of competition kernel}

\item{ac.C}{parameter influencing shape of the function (has to be positive)}

\item{ac.v}{parameter influencing shape of the function (has to be positive)}

\item{...}{Any additional parameters}
}
\description{
It calculates asymmetric competition coefficients
}
\details{
Depending on value of ac.type the convex-concave function from Kisdi (1999)
or  smooth function suggested by Nattrass et al (2012) are used.

For formulas and meaning of parameters see the \code{vignette("competition")}
}
\references{
Kisdi, E. (1999) Evolutionary Branching under Asymmetric
Competition
\emph{Journal of Theoretical Biology} \bold{197}(2): 149-162.
\doi{10.1006/jtbi.1998.0864}

Nattrass, S., Baigent, S., & Murrell, D. J. (2012) Quantifying
the Likelihood of Co-existence for Communities with Asymmetric Competition.
\emph{Bulletin of Mathematical Biology}, \bold{74}(10): 2315–2338.
\doi{10.1007/s11538-012-9755-8}
}
\seealso{
\code{\link{competition.kernel}}
}
