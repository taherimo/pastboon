\name{count_pairwise_transitions}
\alias{count_pairwise_transitions}
\title{
Calculation of pairwise transitions between a given set of states
}
\description{
For a given set of states, the frequencies of transition between each pair of states aftrer some time-steps are calculated.
}
\usage{
count_pairwise_transitions(net, method=c("SDDS","BNp","PEW"), params,
                           states, steps=1, repeats=1000,
                           asynchronous=TRUE, update_prob=NULL)
}
\arguments{
\item{net}{A network structure of class \code{BooleanNetwork} from \pkg{BoolNet}.}

\item{method}{The parameterization method to be used. If "SDDS" is chosen, stochastic discrete dynamical system proposed by Murrugarra et al. will be used. If "BNp" is selected, Boolean network with perturbation initially presented by Golinelli and Derrida (and later by Shmulevich et al. and Trairatphisan et al.) will be utilized. If "PEW" is chosen, Boolean network with probabilistic edge weights proposed by Deritei et al. will be employed. For each of these cases a different pattern for the \code{params} argument is expected.}

\item{params}{The parameter values depending to the value of \code{method}. If \code{method= "SDDS"}, the \code{params} argument is expected to be a list of four equal-size vectors \code{p00}, \code{p01}, \code{p10}, and \code{p11}. If \code{method= "BNp"}, the \code{params} argument is expected to be vector of probabilities. Finally, if \code{params= "PEW"}, the \code{params} argument is expected to be a a list of two equall-size vectors \code{p_on} and \code{p_off}. The length of all probability vectors should be equal to the number of network nodes.}

\iten{states}{A table that each row is binary form of a state. The number of transitions between these states are counted.}

\item{steps}{The number of time steps which should be a non-negatve integer.}

\item{repeats}{The number of time repeats which should be a positive integer.}

\item{asynchronous}{If \code{TRUE}, the asynchronous update scheme is used where a single node is updated at each time-step. In this case \code{update_prob} indicates update probabuilities. If \code{FALSE}, synchronous update is utilized.}

\item{update_prob}{Probability of updating each variable in each time step is the case of \code{asynchronous= TRUE}. If it is set to \code{NULL} (default value) then nodes are randomly updated based on uniform distribution. If \code{asynchronous= FALSE} this argument is ignored.}
}

\details{
Counting the number of transitions between each pair of states is an indicator of reachability of the states.
}

\value{
A matrix that each of its elements i,j shows the number of transitions from the ith state to the jth state.
}

\examples{

# load the example network
data(lac_operon_net)

}
