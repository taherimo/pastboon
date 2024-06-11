\name{get_reached_states}
\alias{get_reached_states}
\title{
Obtain the reached states
}
\description{
Obtains The reached states after simulating a Boolean network for some time-steps.
}
\usage{
get_reached_states(net, method = c("SDDS","BNp","PEW"),
                   params, steps, repeats, initial_states = NULL,
                   asynchronous = TRUE, update_prob = NULL)
}
\arguments{
\item{net}{A network structure of the class \code{BooleanNetwork} from the \pkg{BoolNet} package.}

\item{method}{The parameterization method to be used. If "SDDS", stochastic discrete dynamical system proposed by Murrugarra et al. will be used. If "BNp", Boolean network with perturbation initially presented by Golinelli and Derrida (and later by Shmulevich et al. and Trairatphisan et al.) will be utilized. If "PEW", Boolean network with probabilistic edge weights proposed by Deritei et al. will be employed. For each of these cases a different pattern for the \code{params} argument is expected.}

\item{params}{The parameter values depending on \code{method}. If \code{method = "SDDS"}, it is expected to be a list of four equal-size vectors of probabilities \code{p00}, \code{p01}, \code{p10}, and \code{p11} where all the vector lengths should be equal to the number of network nodes. If \code{method = "BNp"}, it should be a single vector of probabilities with the length equal to the number of network nodes. Finally, if \code{params = "PEW"}, it is required be a list of two equall-size vectors of probabilities \code{p_on} and \code{p_off} with the length equal to the number of network nodes.}

\item{steps}{The number of time-steps, which should be a non-negatve integer.}

\item{repeats}{The number of repeats, which should be a positive integer.}

\item{initial_states}{The list of initial states as a matrix  where the rows are binary form of the states (or single binary vector for a single initial state). he number of matrix columns (the length of the vector) should be equal to the number of network nodes. The  If \code{NULL} the initial state(s) are chosen randomly based on the uniform distrivution.
}

\item{asynchronous}{If \code{TRUE}, the asynchronous update scheme is used where a single node is updated at each time-step. In this case \code{update_prob} indicates update probabuilities. If \code{FALSE}, synchronous update is utilized.}

\item{update_prob}{Probability of updating each variable (node) in each time-step in the case of \code{asynchronous = TRUE}. If \code{NULL} (default) then nodes are randomly updated based on uniform distribution. If \code{asynchronous = FALSE} this argument is ignored.}
}

\details{
This function returns the reached states (the states in the last time-step) after simulating a network for \code{steps} time-steps. If multiple initial states are provided by user, then the \code{repeats} argument is ignored.
}

\value{
The set of reached states.
}

\examples{

# load the example network
data(lac_operon_net)

}

