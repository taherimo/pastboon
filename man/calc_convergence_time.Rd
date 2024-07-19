\name{calc_convergence_time}
\alias{calc_convergence_time}
\title{
Calculate convergence time-step for node activities
}
\description{
Given node activity trajectories for a set of variables \code{node_act}, this function calculates the time-step from which the changes in all the curves are below \code{threshold} for \code{window_size} consecutive time-steps.
}
\usage{
calc_convergence_time(node_act, threshold, window_size = 1)
}

\arguments{
\item{node_act}{A matrix describing node activities over consecutive time-steps, where rows represent time-steps and columns represent nodes. It is the output of \code{\link{calc_node_activities}}.}

\item{threshold}{A value determining the maximum allowable change in node activities to decide if they have converged.}

\item{window_size}{The number of consecutive time-steps for which the node activity curves must remain stable (i.e., changes below \code{threshold}) to be considered converged. The default is 1.}
}
\details{
The function checks if the changes in all node activity curves are less than \code{threshold} for \code{window_size} consecutive time-steps. If this condition is met, the node activity curves are considered to have converged to their stable values, and the convergence time-step (the starting point of the window) is returned. Since node activities represent marginal probabilities of the nodes being active at each time-step, convergence indicates that the steady-state distribution of the corresponding Boolean network has been reached, meaning that the probability of being in each state of the network no longer changes significantly.
}

\value{
The time-step at which convergence occurs. If no convergence is detected, \code{NA} is returned.
}

\examples{

# Load the example network
data(lac_operon_net)

# Define parameters for the SDDS method
props <- rep(0.95, length(lac_operon_net$genes))
params <- list(p00 = props, p01 = props, p10 = props, p11 = props)

# Get node activities after simulation using the SDDS method
node_act <- calc_node_activities(lac_operon_net, method = "SDDS", params = params,
    steps = 100, repeats = 10000)

# Calculate the convergence time
convergence_time <- calc_convergence_time(node_act, threshold = 0.01)


# Print the convergence time
print(convergence_time)
}
