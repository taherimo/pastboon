\name{count_pairwise_trans}
\alias{count_pairwise_trans}
\title{
Count pairwise transitions between a given set of states
}
\description{
Counts the frequencies of transitions between each pair of states from a given set of states.
}
\usage{
count_pairwise_trans(net, method = c("BNp", "SDDS", "PEW"), params, states,
                     steps = 1, repeats = 1000, asynchronous = TRUE,
                     update_prob = NULL)
}
\arguments{
\item{net}{A network structure of the class \code{BooleanNetwork} from the \code{BoolNet} package.}

\item{method}{The parameterization method to be used. Options are:

  \itemize{

    \item \code{"BNp"}: Boolean network with perturbations.

    \item \code{"SDDS"}: Stochastic discrete dynamical systems.

    \item \code{"PEW"}: Boolean network with probabilistic edge weights.

  }

  Each method requires a different format for the \code{params} argument.
}

\item{params}{The parameter values depending on \code{method}:
\itemize{
    \item For \code{method = "BNp"}, a single vector of probabilities, equal in length to the number of network nodes.

    \item For \code{method = "SDDS"}, a list of four equal-length vectors of probabilities: \code{p00}, \code{p01}, \code{p10}, and \code{p11}, each equal in length to the number of network nodes.

    \item For \code{method = "PEW"}, a list of two equal-length vectors of probabilities: \code{p_on} and \code{p_off}, each as long as the number of network edges, ordered according to \code{\link{extract_edges}}.
}
}

\item{states}{The network states among which pairwise transitions are to be counted. This should be a matrix (where the rows represent the binary form of the states) or a vector (for the binary form of a single state). The number of matrix columns (or the length of the vector) should match the number of network nodes.}

\item{steps}{The number of time-steps, which should be a non-negative integer.}

\item{repeats}{The number of repeats, which should be a positive integer.}

\item{asynchronous}{If \code{TRUE}, the asynchronous update scheme is used, where a single node is updated at each time-step. In this case, \code{update_prob} indicates update probabilities. If \code{FALSE}, the synchronous update scheme is utilized.}

\item{update_prob}{The probability of updating each variable (node) in each time-step when \code{asynchronous = TRUE}. If \code{NULL} (default), nodes are updated randomly based on a uniform distribution. If \code{asynchronous = FALSE}, this argument is ignored.}
}

\details{
Counting the number of transitions between each pair of states reveals the reachability of one state from another. This function performs simulations by starting from each state in \code{states} for \code{steps} time-steps and \code{repeats} iterations, and counts the number of transitions to other states in \code{states}.
}

\value{
A matrix where each element (i, j) represents the number of transitions from the \emph{i}th state to the \emph{j}th state across \code{steps} time-steps and \code{repeats} iterations.
}

\references{
Golinelli, O., & Derrida, B. (1989). Barrier heights in the Kauffman model. \emph{Journal De Physique}, 50(13), 1587-1601.

Shmulevich, I., Dougherty, E. R., & Zhang, W. (2002). Gene perturbation and intervention in probabilistic Boolean networks. \emph{Bioinformatics}, 18(10), 1319-1331.

Trairatphisan, P., Mizera, A., Pang, J., Tantar, A. A., Schneider, J., & Sauter, T. (2013). Recent development and biomedical applications of probabilistic Boolean networks. \emph{Cell communication and signaling}, 11, 1-25.

Murrugarra, D., Veliz-Cuba, A., Aguilar, B., Arat, S., & Laubenbacher, R. (2012). Modeling stochasticity and variability in gene regulatory networks. \emph{EURASIP Journal on Bioinformatics and Systems Biology}, 2012, 1-11.

Deritei, D., Kunšič, N., & Csermely, P. (2022). Probabilistic edge weights fine-tune Boolean network dynamics. \emph{PLoS Computational Biology}, 18(10), e1010536.
}

\examples{

# >>>>>>>>>>>>>>>>  Load network and generate random initial states  <<<<<<<<<<<<<<<<<

# Load the example network
data(lac_operon_net)

# Generate some random states
states <- matrix(sample(c(0, 1), 10 * length(lac_operon_net$genes),
    replace = TRUE), nrow = 10, ncol = length(lac_operon_net$genes))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Method: BNp  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Define the parameters for the BNp method
params <- rep(0.05, length(lac_operon_net$genes))

# Obtain frequency of pairwise transitions
pairwise_trans <- count_pairwise_trans(lac_operon_net, states = states,
    method = "BNp", params = params, steps = 100, repeats = 10)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Method: SDDS  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Define the parameters for the SDDS method
props <- rep(0.95, length(lac_operon_net$genes))
params <- list(p00 = props, p01 = props, p10 = props, p11 = props)

# Obtain frequency of pairwise transitions
pairwise_trans <- count_pairwise_trans(lac_operon_net, states = states,
    method = "SDDS", params = params, steps = 100, repeats = 10)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Method: PEW  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Extract edges from the network
edges <- extract_edges(lac_operon_net)

# Define the parameters for the PEW method
p_on <- runif(nrow(edges))
p_off <- runif(nrow(edges))
params <- list(p_on = p_on, p_off = p_off)

# Obtain frequency of pairwise transitions
pairwise_trans <- count_pairwise_trans(lac_operon_net, states = states,
    method = "PEW", params = params, steps = 100, repeats = 10)
}
