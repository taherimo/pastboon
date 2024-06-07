\name{extract_edges}
\alias{extract_edges}
\title{
Extracts edges from a Boolean network
}
\description{
The list of directed edges (links) are extracted from a given Boolean network.
}
\usage{
extract_edges(network, node_names= TRUE)
}
\arguments{
\item{net}{A network structure of class \code{BooleanNetwork} from \pkg{BoolNet}.}

\item{node_names}{If TRUE (default), the edges are returned based on node names, otherwise the edges are retuend based on node indices.}
}

\details{
Since the topology of each Boolean network is a directed graph, this function extracts the list of directed edges from a given Boolean network.
}

\value{
A table is returned where its rows correspond to the directed edges of the network and its two columns idicate source and destination of each edge.
}

\examples{

# load the example network
data(lac_operon_net)

edges <- extract_edges(lac_operon_net)

}
