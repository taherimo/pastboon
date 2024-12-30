\name{myeloid_differentiation_net}

\alias{myeloid_differentiation_net}

\docType{data}

\title{
The myeloid differentiation Boolean network
}

\description{
The myeloid differentiation Boolean network proposed by Krumsiek et al.
}

\usage{data(myeloid_differentiation_net)}
\details{
The data consists of an object \code{myeloid_differentiation_net} of the class \code{BooleanNetwork} (from the \code{BoolNet} package), describing the myeloid differentiation gene regulatory network with 11 genes. More specifically, this network represents differentiation of common myeloid progenitors to megakaryocytes, erythrocytes, granulocytes and monocytes. The state space of the model is a hierarchical, acyclic graph, showing the principles of myeloid differentiation. A good agreement between the steady states of the model and microarray expression profiles of two different studies has been observed.
}

\references{
Krumsiek, J., Marr, C., Schroeder, T., & Theis, F. J. (2011). Hierarchical differentiation of myeloid progenitors is encoded in the transcription factor network. \emph{PloS one}, 6(8), e22649.
}

\examples{
# load the network
data(myeloid_differentiation_net)

# the network is stored in a variable called 'myeloid_differentiation_net'
print(myeloid_differentiation_net)
}
