\name{read.OCG}
\alias{read.OCG}

\title{
Read an OCG Partition File into R
}
\description{
This function reads in an OCG partition file and converts it into an OCG object for use in R.
}
\usage{
read.OCG(file, elfile = NULL, verbose = FALSE, keep.out = FALSE)
}

\arguments{
  \item{file}{
A character string naming the OCG partition file.
}
  \item{elfile}{
A character string naming the file containing the network that the OCG partition is based upon.
}
  \item{verbose}{
Logical, whether to print progress to the screen. Defaults to FALSE.
}
  \item{keep.out}{
Logical, whether to keep the intermediate files written when reading in the OCG partition. Defaults to FALSE.
}
}

\value{
An object of class \code{OCG}, which is a list containing the following elements:
\item{numbers }{An integer vector with the number of edges, nodes, and communities.}
\item{modularity }{An integer number specifying the modularity of the network.}
\item{Q }{A real number specifying the value of Q generated by the OCG algorithm.}
\item{nodeclusters }{A data frame consisting of 2 columns; the first contains node names, and the second contains single community IDs for each node. All communities and their nodes are represented, but not necessarily all nodes.}
\item{numclusters }{A named integer vector. Names are node names and integer values are the number of communities to which each node belongs.}
\item{igraph }{An object of class \code{\link{igraph}}. The network is represented here as an \code{igraph} object.}
\item{edgelist }{A character matrix with 2 columns containing the nodes that interact with each other.}
\item{clustsizes }{A named integer vector. Names are community IDs and integer values indicate the number of nodes that belong in each community.}
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\examples{
## Read an OCG partition file into R.
\dontrun{oc <- read.OCG(file = "OCG_partition.txt", elfile = "network.txt")}
}

