\name{get_reached_states}
\alias{get_reached_states}
\title{
Obtain the reached states
}
\description{
Obtains the reached states after simulating a Boolean network for some time-steps.
}
\usage{
get_reached_states(net, method = c("SDDS","BNp","PEW"), params, steps,
                   repeats = NULL, initial_states = NULL, asynchronous = TRUE,
                   update_prob = NULL)
}
\arguments{
\item{net}{A network structure of the class \code{BooleanNetwork} from the \pkg{BoolNet} package.}

\item{method}{The parameterization method to be used. For "SDDS", "BNp", amd
"PEW" the methods \emph{stochastic discrete dynamical systems}, \emph{Boolean network with perturbations}, and \emph{Boolean network with probabilistic edge weights} will be employed, respectively. For each of these cases a different pattern for the \code{params} argument is expected. For more information, see \code{\link{calc_node_activities}}.
}

\item{params}{The parameter values depending on \code{method}. If \code{method = "SDDS"}, it is expected to be a list of four equal-size vectors of probabilities \code{p00}, \code{p01}, \code{p10}, and \code{p11} where all the vector lengths should be equal to the number of network nodes. If \code{method = "BNp"}, it should be a single vector of probabilities with the length equal to the number of network nodes. Finally, if \code{params = "PEW"}, it is required be a list of two equall-size vectors of probabilities \code{p_on} and \code{p_off} with the length equal to the number of network edges based on the order of the edges obtained from \code{\link{extract_edges}}.}

\item{steps}{The number of time-steps, which should be a non-negatve integer.}

\item{repeats}{The number of repeats, which should be a positive integer. If two or more initial states are provided via \code{initial_states}, then this argument is ignored. If \code{NULL} (default), then \code{initial_states} should not be \code{NULL}.}

\item{initial_states}{The set of initial states as a matrix  where each row corresponds to the binary form of a state (or single binary vector for a single initial state). The number of matrix columns (the length of the vector) should be equal to the number of network nodes. The order of the nodes in the columns is considerd same as \code{net$genes}. If \code{NULL} (default), the initial state(s) are chosen randomly for \code{repeats} number of times based on the uniform distribution, thus \code{repeats} should not be \code{NULL} ate the same time.
}

\item{asynchronous}{If \code{TRUE}, the asynchronous update scheme is used where a single node is updated at each time-step. In this case \code{update_prob} indicates update probabuilities. If \code{FALSE}, synchronous update is utilized.}

\item{update_prob}{Probability of updating each variable (node) in each time-step in the case of \code{asynchronous = TRUE}. If \code{NULL} (default) then nodes are randomly updated based on uniform distribution. If \code{asynchronous = FALSE} this argument is ignored.}
}

\details{
This function returns the reached states (the states in the last time-step) after simulating a network for \code{steps} time-steps and repeating it for \code{repeats} number of times. If \code{initial_states = NULL}, then the intial states are chosen randomly based on uniform distribution for \code{repeats} number of times and as a result, \code{repeats} number of reached states are returned. If two or more initial states are provided by user, then the \code{repeats} argument is ignored, and one reached state is returned for each initial state. If \code{repeats = NULL}, then the number of the returned reached states is equal to the number of inirial states. i.e. for each initial state a reached state is returned. Please note that \code{repeats} and \code{initial_states} should not be \code{NULL} at the same time.
}

\value{
A matrix where each row is the binary form of a reached state and each column corresponds to a network node. The order of the nodes in the columns is same as \code{net$genes}.
}

\references{
Murrugarra, D., Veliz-Cuba, A., Aguilar, B., Arat, S., & Laubenbacher, R. (2012). Modeling stochasticity and variability in gene regulatory networks. \emph{EURASIP Journal on Bioinformatics and Systems Biology}, 2012, 1-11.

Golinelli, O., & Derrida, B. (1989). Barrier heights in the Kauffman model. \emph{Journal De Physique}, 50(13), 1587-1601.

Shmulevich, I., Dougherty, E. R., & Zhang, W. (2002). Gene perturbation and intervention in probabilistic Boolean networks. \emph{Bioinformatics}, 18(10), 1319-1331.

Trairatphisan, P., Mizera, A., Pang, J., Tantar, A. A., Schneider, J., & Sauter, T. (2013). Recent development and biomedical applications of probabilistic Boolean networks. \emph{Cell communication and signaling}, 11, 1-25.

Deritei, D., Kunšič, N., & Csermely, P. (2022). Probabilistic edge weights fine-tune Boolean network dynamics. \emph{PLoS Computational Biology}, 18(10), e1010536.
}

\examples{

# >>>>>>>>>>>>>>>>  Load network and generate random initial states  <<<<<<<<<<<<<<<<<

# Load the example network
data(lac_operon_net)

# Generate a single initial state
initial_state <- sample(c(0, 1), length(lac_operon_net$genes), replace = TRUE)

# Generate multiple (10) initial states
initial_states <- matrix(sample(c(0, 1), 10 * length(lac_operon_net$genes),
    replace = TRUE), nrow = 10, ncol = length(lac_operon_net$genes))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Method: SDDS  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Define the parameters for the SDDS method
props <- rep(0.95, length(lac_operon_net$genes))
params <- list(p00 = props, p01 = props, p10 = props, p11 = props)

# No initial states are provided
reached_states <- get_reached_states(lac_operon_net, method = "SDDS", params = params,
    steps = 100, repeats = 10)

# A single initial state is provided
reached_states <- get_reached_states(lac_operon_net, method = "SDDS", params = params,
    steps = 100, initial_states = initial_state, repeats = 10)

# Multiple initial states are provided
reached_states <- get_reached_states(lac_operon_net, method = "SDDS", params = params,
    steps = 100, initial_states = initial_states)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Method: BNp  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Define the parameters for the BNp method
params <- rep(0.05, length(lac_operon_net$genes))

# No initial states are provided
reached_states <- get_reached_states(lac_operon_net, method = "BNp", params = params,
    steps = 100, repeats = 10)

# A single initial state is provided
reached_states <- get_reached_states(lac_operon_net, method = "BNp", params = params,
    steps = 100, initial_states = initial_state, repeats = 10)

# Multiple initial states are provided
reached_states <- get_reached_states(lac_operon_net, method = "BNp", params = params,
    steps = 100, initial_states = initial_states)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Method: PEW  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Extract edges from the network
edges <- extract_edges(lac_operon_net)

# Define the parameters for the PEW method
p_on <- runif(nrow(edges))
p_off <- runif(nrow(edges))
params <- list(p_on = p_on, p_off = p_off)

# No initial states are provided
reached_states <- get_reached_states(lac_operon_net, method = "PEW", params = params,
    steps = 100, repeats = 10)

# A single initial state is provided
reached_states <- get_reached_states(lac_operon_net, method = "PEW", params = params,
    steps = 100, initial_states = initial_state, repeats = 10)

# Multiple initial states are provided
reached_states <- get_reached_states(lac_operon_net, method = "PEW", params = params,
    steps = 100, initial_states = initial_states)
}

