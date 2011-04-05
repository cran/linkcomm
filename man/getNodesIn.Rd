\name{getNodesIn}
\alias{getNodesIn}

\title{
Extract Nodes from Communities
}
\description{
This function returns node names that belong to sets of communities.
}
\usage{
getNodesIn(x, clusterids = 1)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{clusterids}{
An integer vector of community IDs. Defaults to community 1.
}
}

\value{
A character vector of node names.
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

## Get nodes from communities 1 and 3.
getNodesIn(lc, clusterids = c(1,3))
}

