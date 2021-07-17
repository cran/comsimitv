\name{tolerance}
\alias{tolerance}
\title{Habitat suitability (tolerance) functions}

\arguments{
\item{trait.values}{Values of trait related to habitat filtering}
\item{env}{Vector of environmental conditions in the local communities}
\item{...}{Additional parameters}
}
\description{
User defined functions for habitat suitability
}
\value{
A matrix of habitat suitabilities, communities in rows, individuals in columns
}
\details{
User can defined any specific function of habitat suitability, depending on environmental conditions and trait value related to habitat filtering.
Vectors of these trait values and environmental conditions have to be
the first and second parameter of the function, and any further parameters are allowed.
The output has to be a matrix of habitat suitabilities, communities in rows, species/individuals in columns.

Tolerance functions available in the package:

\code{\link{Gaussian.tolerance}}
}
