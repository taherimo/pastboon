\name{extract_edges}
\alias{extract_edges}
\title{
Extract edges from a Boolean network
}
\description{
Extracts the list of directed edges (links) from a given Boolean network.
}
\usage{
extract_edges(net, node_names = TRUE)
}
\arguments{
\item{net}{A network structure of the class \code{BooleanNetwork} from the \code{BoolNet} package.}

\item{node_names}{If TRUE (default), the edges are returned by node names; otherwise, the edges are returned by node indices.}
}

\details{
Since Boolean networks have a directed graph topology, this function extracts the list of directed edges from a given Boolean network.
}

\value{
A data frame where each row corresponds to a directed edge of the network and the two columns indicate the source and destination of each edge.
}

\examples{

# Load the example network
data(lac_operon_net)

# Extract edges from the network
edges <- extract_edges(lac_operon_net)
}
