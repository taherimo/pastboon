\name{get_convergence_time}
\alias{get_convergence_time}
\title{
Calculate the time that a set of node activity curves are converged
}
\description{
Given node activity trajectories for a set of variables \code{node_activities}, this function calculates the time-step that the changes in all the curves are below \cite{threshold} in the next \code{window_size} consecutive time-steps.
}
\usage{
get_convergence_time(node_activities,
                     threshold,
                     window_size=1)
}

\arguments{
  \item{node_activiries}{
A table that describes node activities in different time-steps where its rows and columns correspond to the time-steps and the nodes, respectively (the output of \code{\link{loadNetwork}}).
}
\item{threshold}{
A value that determines the maximum change in node activiies to decide if they are concerged.
}
\item{window_size}{
The number of the consecutive time-steps that the curves are required to be stable (changes are below \cite{threshold}). The default is 1.
}
}
\details{

}
\examples{
# load the example network
data(lac_operon_net)

}
