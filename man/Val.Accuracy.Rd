% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Val.Accuracy.R
\name{Val.Accuracy}
\alias{Val.Accuracy}
\title{Network Inference Accuracy Validation}
\usage{
Val.Accuracy(Network1, Network2)
}
\arguments{
\item{Network1}{igraph object. First network in comparison.}

\item{Network2}{igraph object. Second network in comparison.}
}
\value{
Percentage of shared links between the two networks out of all possible links.
}
\description{
Calculate accuracy of network inference by tallying between two networks the number of similar links by their sign (i.e., absent, positive, and negative links) and divide the number of shared links by the number of all possible links. Networks must be of the same size and contain the same node sets.
}
\examples{
data("Network_igraph")
data("Inferred_igraph")
Val.Accuracy(Network1 = Network_igraph, Network2 = Inferred_igraph)

}
