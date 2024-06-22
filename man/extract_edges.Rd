\name{extract_edges}
\alias{extract_edges}
\title{
Extract edges from a Boolean network
}
\description{
Extracts the list of directed edges (links) from a given Boolean network.
}
\usage{
extract_edges(network, node_names = TRUE)
}
\arguments{
\item{net}{A network structure of the class \code{BooleanNetwork} from the \pkg{BoolNet} package.}

\item{node_names}{If TRUE (default), the edges are returned by node names, otherwise the edges are retuend by node indices.}
}

\details{
Since the Boolean network have a directed graph topology, this function extracts the list of directed edges from a given Boolean network.
}

\value{
A data frame is returned where its rows correspond to the directed edges of the network and its two columns idicate source and destination of each edge.
}

\examples{

# Load the example network
data(lac_operon_net)

# Extract edges from the network
edges <- extract_edges(lac_operon_net)
}
