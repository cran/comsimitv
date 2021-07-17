\name{fITV}
\alias{fITV}
\title{Intraspecific Trait Variation}

\arguments{
\item{seeds}{Matrix of produced seeds (with mother's trait values) as produced by
\code{\link{SeedProduction}} function}

\item{...}{Other parameters of the function}

}
\value{
The same type as \code{seeds} parameter, i.e. a matrix where each seed is represented by one row, and seeds's attributes (location, species, trait values) are in the columns
}
\description{
User defined function for Intraspecific Trait Variation
}
\details{
User can defined any specific function for ITV, e.g. random variation around mothers value, or maternal effect.

The first parameter has to be matrix of produced seeds, in the form as it created by \code{\link{SeedProduction}} function, and the results has to be in the same matrix form with updated trait values.


ITV functions available in the package:

\code{\link{randomITV}}
}



