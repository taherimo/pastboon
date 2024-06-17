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

\item{method}{The parameterization method to be used. For "SDDS", "BNp", amd
"PEW" the methods \emph{stochastic discrete dynamical systems}, \emph{Boolean network with perturbations}, and \emph{Boolean network with probabilistic edge weights} will be employed, respectively. For each of these cases a different pattern for the \code{params} argument is expected. For more information, see \code{\link{calc_node_activities}}.
}

\item{params}{The parameter values depending on \code{method}. If \code{method = "SDDS"}, it is expected to be a list of four equal-size vectors of probabilities \code{p00}, \code{p01}, \code{p10}, and \code{p11} where all the vector lengths should be equal to the number of network nodes. If \code{method = "BNp"}, it should be a single vector of probabilities with the length equal to the number of network nodes. Finally, if \code{params = "PEW"}, it is required be a list of two equall-size vectors of probabilities \code{p_on} and \code{p_off} with the length equal to the number of network edges based on the order of the edges obtained from \code{\link{extract_edges}}.}

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

\references{
Murrugarra, D., Veliz-Cuba, A., Aguilar, B., Arat, S., & Laubenbacher, R. (2012). Modeling stochasticity and variability in gene regulatory networks. \emph{EURASIP Journal on Bioinformatics and Systems Biology}, 2012, 1-11.
\cr
\cr
Deritei, D., Kunšič, N., & Csermely, P. (2022). Probabilistic edge weights fine-tune Boolean network dynamics. \emph{PLoS Computational Biology}, 18(10), e1010536.
\cr
\cr
Golinelli, O., & Derrida, B. (1989). Barrier heights in the Kauffman model. \emph{Journal De Physique}, 50(13), 1587-1601.

Shmulevich, I., Dougherty, E. R., & Zhang, W. (2002). Gene perturbation and intervention in probabilistic Boolean networks. \emph{Bioinformatics}, 18(10), 1319-1331.

Trairatphisan, P., Mizera, A., Pang, J., Tantar, A. A., Schneider, J., & Sauter, T. (2013). Recent development and biomedical applications of probabilistic Boolean networks. \emph{Cell communication and signaling}, 11, 1-25.
}

\examples{

# load the example network
data(lac_operon_net)

# /\/\/\/\ Example 1: The SDDS Method /\/\/\/\

# define the parameter for the SDDS method
props <- rep(0.95, length(lac_operon_net$genes))
params <- list(p00= props, p01= props, p10= props, p11= props)

reached_states <- count_pairwise_trans(net, method = "SDDS", params = params, steps = 100, repeats = 10)


# /\/\/\/\ Example 2: The BNp Method /\/\/\/\


# define the parameters for the BNp method
params <- rep(0.05, length(lac_operon_net$genes))


# /\/\/\/\ Example 3: The PEW Method /\/\/\/\


# extract edges from the network
edges <- extract_edges(lac_operon_net)

# define the parameters for the PEW method
p_on <- runif(nrow(edges))
p_off <- runif(nrow(edges))
params <- list(p_on=p_on, p_off=p_off)

}



