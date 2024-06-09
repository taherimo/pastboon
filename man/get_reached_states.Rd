\name{get_reached_states}
\alias{get_reached_states}
\title{
Returns the reached states
}
\description{
The reached states are obtained after simulating a Boolean network for some time-steps.
}
\usage{
get_reached_states(net, method=c("SDDS","BNp","PEW"),
                   params, steps, repeats, initial_states= NULL,
                   update_prob= NULL,asynchronous= TRUE)
}
\arguments{
\item{net}{A network structure of class \code{BooleanNetwork} from \pkg{BoolNet}.}

\item{method}{The parameterization method to be used. If "SDDS" is chosen, stochastic discrete dynamical system proposed by Murrugarra et al. will be used. If "BNp" is selected, Boolean network with perturbation initially presented by Golinelli and Derrida (and later by Shmulevich et al. and Trairatphisan et al.) will be utilized. If "PEW" is chosen, Boolean network with probabilistic edge weights proposed by Deritei et al. will be employed. For each of these cases a different pattern for the \code{params} argument is expected.}

\item{params}{The parameter values depending to the value of \code{method}. If \code{method= "SDDS"}, the \code{params} argument is expected to be a list of four equal-size vectors \code{p00}, \code{p01}, \code{p10}, and \code{p11}. If \code{method= "BNp"}, the \code{params} argument is expected to be vector of probabilities. Finally, if \code{params= "PEW"}, the \code{params} argument is expected to be a a list of two equall-size vectors \code{p_on} and \code{p_off}. The length of all probability vectors should be equal to the number of network nodes.}

\item{steps}{The number of time-steps which should be a non-negatve integer.}

\item{repeats}{The number of time repeats which should be a positive integer.}

\item{asynchronous}{If \code{TRUE}, the asynchronous update scheme is used where a single node is updated at each time-step. In this case \code{update_prob} indicates update probabuilities. If \code{FALSE}, synchronous update is utilized.}
}
