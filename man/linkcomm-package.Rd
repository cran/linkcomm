\name{linkcomm-package}
\alias{linkcomm-package}
\alias{linkcomm}
\docType{package}
\title{
The \code{linkcomm} package}
\description{
\code{linkcomm} provides tools for the generation, visualization, and analysis of link communities in networks of arbitrary size and type.
}

\details{
Link communities reveal the nested and overlapping structure in networks, and uncover the key nodes that form connections to multiple communities. \code{linkcomm} provides tools for generating, visualizing, and analysing link communities in networks of arbitrary size and type.

For a more detailed overview of how to use the package:\cr \cr
\code{vignette(topic = "linkcomm", package = "linkcomm")}\cr \cr
To run an interactive demonstration of \code{linkcomm} within \R:\cr \cr
\code{demo(topic = "linkcomm", package = "linkcomm")}
}

\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}
\references{
Ahn, Y.Y., Bagrow, J.P., and Lehmann, S. (2010). Link communities reveal multiscale complexity in networks. Nature 466, 761-764.\cr
Spencer, R. (2010). \url{http://scaledinnovation.com/analytics/communities/comlinks.html}
}

\seealso{
\code{\link{getLinkCommunities}}, \code{\link{plot.linkcomm}}, \code{\link{pp_rnapol}}, \code{\link{lesmiserables}}, \code{\link{karate}}, \code{\link{weighted}}, \code{\link{igraph}}, \code{\link{RColorBrewer}}, \code{\link{grid}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Plot a graph layout of the link communities.
plot(lc, type = "graph")

## Use a Spencer circle layout.
plot(lc, type = "graph", layout = "spencer.circle")

## Calculate a community-based measure of node centrality.
getCommunityCentrality(lc)

## Find nested communities.
getAllNestedComm(lc)

## Uncover the relatedness between communities.
getClusterRelatedness(lc)
}

