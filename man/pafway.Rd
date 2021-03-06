% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mainpafway.R
\name{pafway}
\alias{pafway}
\title{Find pairwise-associations between annotations in a network without edge weights.}
\usage{
pafway(GO, edges, GOtypes, exact = TRUE, adjustByEdgeCount = FALSE)
}
\arguments{
\item{GO}{A vector of Strings, equal to the length of all the nodes.  The names of the vector should be the names of the nodes.  The values should either be the functional annotations, or a concatenation of the functional annotations, separated by a "_" symbol.}

\item{edges}{A matrix of Strings, with at least two columns.  Each row will represent an edge, linking the node in the first column to the node in the second column.  Please make sure the node names are the same as those in "GO"}

\item{GOtypes}{This is a vector that contains the functional annotations or GO terms that are of interest}

\item{exact}{A boolean.  If it is true, it will look for an exact match between the term in GOtypes and the vector GO.  Otherwise, it will look for substrings.}

\item{adjustByEdgeCount}{A boolean.  If true, then the probability of observing a functional annotation will be calculated in terms of the number of nodes, but if it is false, it is calculated in terms of the number of edges that contain that node.}
}
\value{
A matrix that has the same number of rows and columns as length(GOtypes).  This will contain p-values.
}
\description{
Find pairwise-associations between annotations in a network without edge weights.
}
\examples{
nodes=paste("node", c(1:10))
set.seed(123)
randomGO=c("A", "B", "C")[sample(c(1:3), 10, replace=TRUE)]
names(randomGO)=nodes
edgesRandom=sapply(c(1:100), function(i){
   nodes[sample(10, 2)]
})
pafway(randomGO, t(edgesRandom), unique(randomGO))
}
