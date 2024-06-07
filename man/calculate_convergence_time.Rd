\name{calculate_convergence_time}
\alias{calculate_convergence_time}
\title{
Calculates converegence time for node activity curve(s)
}
\description{
Given node activity curves for a set of variables \code{node_activities}, this function calculates the time-step that the changes in all the curves are below \cite{threshold} in the next \code{window_size} consecutive time-steps.
}
\usage{
get_convergence_time(node_activities, threshold, window_size=1)
}

\arguments{
\item{node_activities}{A table that describes node activities in different time-steps where its rows and columns correspond to the time-steps and nodes, respectively. It is the output of \code{\link{calculate_node_activities}}.
}
\item{threshold}{A value that determines the maximum change in node activiies to decide if they are concerged.}

\item{window_size}{The number of the consecutive time-steps that the curves are required to be stable (changes are below \cite{threshold}). The default is 1.}
}
\details{
If the change of all curves are less than \code{threshold} then it is assumed that the the convergene time step is reached.
}
\examples{
# load the example network
data(lac_operon_net)

# define the parameter for the SDDS method
props <- rep(0.95, length(lac_operon_net$genes)
params <- list(p00= props, p01= props, p10= props, p11= props)

# get node activities after simulation using the SDDS method and asynchronous update
node_act <- calculate_node_activities(lac_operon_net, method= "SDDS", params =                                            params, steps= 100, repeats = 10000)

calculate_convergence_time(node_act)

}
