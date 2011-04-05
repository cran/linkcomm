\name{getCommunityMatrix}
\alias{getCommunityMatrix}

\title{
Construct a Community Membership Matrix
}
\description{
This function returns a binary matrix with nodes as rows, communities as columns, and unit entries indicating membership in a community.
}
\usage{
getCommunityMatrix(x, nodes = head(names(x$numclusters), 20))
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{nodes}{
A character vector containing the nodes for the community membership matrix. Defaults to the 20 (or less) nodes that belong to the most communities.
}
}

\value{
A binary matrix with nodes as rows and communities as columns.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{plot.linkcomm}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Get community membership matrix.
getCommunityMatrix(lc)
}

