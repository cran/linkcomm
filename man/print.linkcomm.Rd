\name{print.linkcomm}
\alias{print.linkcomm}

\title{
Print a Summary of a \code{linkcomm} Object
}
\description{
This function prints summary statistics for a \code{linkcomm} object to the screen.
}
\usage{
\method{print}{linkcomm}(x, ...)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{\dots}{
Further arguments passed to or from other methods.
}
}

\value{
Prints summary data to the screen.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Print summary statistics to the screen.
print(lc)
}

