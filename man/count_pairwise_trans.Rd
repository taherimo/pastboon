\name{count_pairwise_trans}
\alias{count_pairwise_trans}
\title{
Count pairwise transitions between a given set of states
}
\description{
Counts the frequencies of transition between each pair of states from a given set of states.
}
\usage{
count_pairwise_trans(net, method = c("SDDS","BNp","PEW"), params,
                           states, steps = 1, repeats = 1000,
                           asynchronous = TRUE, update_prob = NULL)
}
\arguments{
\item{net}{A network structure of the class \code{BooleanNetwork} from the \pkg{BoolNet} package.}

\item{method}{The parameterization method to be used. If "SDDS", stochastic discrete dynamical system proposed by Murrugarra et al. will be used. If "BNp", Boolean network with perturbation initially presented by Golinelli and Derrida (and later by Shmulevich et al. and Trairatphisan et al.) will be utilized. If "PEW", Boolean network with probabilistic edge weights proposed by Deritei et al. will be employed. For each of these cases a different pattern for the \code{params} argument is expected.}

\item{params}{The parameter values depending on \code{method}. If \code{method = "SDDS"}, it is expected to be a list of four equal-size vectors of probabilities \code{p00}, \code{p01}, \code{p10}, and \code{p11} where all the vector lengths should be equal to the number of network nodes. If \code{method = "BNp"}, it should be a single vector of probabilities with the length equal to the number of network nodes. Finally, if \code{params = "PEW"}, it is required be a list of two equall-size vectors of probabilities \code{p_on} and \code{p_off} with the length equal to the number of network nodes.}

\iten{states}{The network states that pairwise transitions among them are going to be counted. Should be a matrix where the rows are binary form of the states (or single binary vector for a single state). The number of matrix columns (the length of the vector) should be equal to the number of network nodes.}

\item{steps}{The number of time-steps, which should be a non-negatve integer.}

\item{repeats}{The number of repeats, which should be a positive integer.}

\item{asynchronous}{If \code{TRUE}, the asynchronous update scheme is used where a single node is updated at each time-step. In this case \code{update_prob} indicates update probabuilities. If \code{FALSE}, synchronous update is utilized.}

\item{update_prob}{Probability of updating each variable (node) in each time-step in the case of \code{asynchronous = TRUE}. If \code{NULL} (default) then nodes are randomly updated based on uniform distribution. If \code{asynchronous = FALSE} this argument is ignored.}
}

\details{
Counting the number of transitions between each pair of states reveals the reachability of one state from another. This function performs simulations by starting from each state in \code{states} for \code{steps} time-steps and \code{repeats} number of times, and counts the number of transitions to other states in \code{states}.
}

\value{
A matrix where each element (i, j) shows the number of transitions from the \emph{i}th state to the \emph{j}th state across \code{steps} time-steps and \code{repeats}.
}

\examples{

# load the example network
data(lac_operon_net)

}
