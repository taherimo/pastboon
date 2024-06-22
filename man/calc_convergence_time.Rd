\name{calc_convergence_time}
\alias{calc_convergence_time}
\title{
Calculate converegence time-step for node activities
}
\description{
Given node activity trajectories for a set of variables \code{node_act}, this function calculates the time-step that the changes in all the curves are below \code{threshold} in the next \code{window_size} consecutive time-steps.
}
\usage{
get_convergence_time(node_act, threshold, window_size = 1)
}

\arguments{
\item{node_act}{A matrix that describes node activities in different time-steps where its rows and columns correspond to the time-steps and nodes, respectively. It is the output of \code{\link{calc_node_activities}}.
}
\item{threshold}{A value that determines the maximum change in node activiies to decide if they are concerged.}

\item{window_size}{The number of the consecutive time-steps that the curves are required to be stable (changes are below \code{threshold}). The default is 1.}
}
\details{
If the change of all curves are less than \code{threshold} then it is assumed that the the convergene time step is reached.
}

\value{
The time-step at which the curves are converged. If no convergence time is detected \code{NA} is returned.
}

\examples{

# load the example network
data(lac_operon_net)

# define the parameter for the SDDS method
props <- rep(0.95, length(lac_operon_net$genes)
params <- list(p00= props, p01= props, p10= props, p11= props)

# get node activities after simulation using the SDDS method and asynchronous update
node_act <- calc_node_activities(lac_operon_net, method= "SDDS", params = params,
    steps= 100, repeats = 10000)

calc_convergence_time(node_act)
}
