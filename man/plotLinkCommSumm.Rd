\name{plotLinkCommSumm}
\alias{plotLinkCommSumm}

\title{
Plot a Summary of the Link Community Algorithm Output
}
\description{
This function is called by \code{plot.linkcomm} to plot a summary of the output of the \code{linkcomm} algorithm.
}
\usage{
plotLinkCommSumm(x, col = TRUE, pal = brewer.pal(9, "Set1"), right = TRUE, 
                 droptrivial = TRUE, ...)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{col}{
Logical, whether to colour link communities in the dendrogram. Defaults to TRUE.
}
  \item{pal}{
A character vector describing a colour palette to be used for colouring the link community dendrogram. Defaults to \code{brewer.pal(9, "Set1")}.
}
  \item{right}{
Logical, whether to orient the dendrogram to the right. Defaults to TRUE.
}
  \item{droptrivial}{
Logical, whether to not colour communities of size 2. Defaults to TRUE.
}
  \item{\dots}{
Additional arguments to be passed to \code{plot}.
}
}
\details{
Here we describe the parameters for plotting link community summaries using:\cr
\code{plot(x, type = "summary")}
}
\value{
A summary plot of the output from the \code{linkcomm} algorithm.
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

## Plot the modularity of the link communities.
plot(lc, type = "summary")
}

