\name{calc_node_activities}
\alias{calc_node_activities}
\title{
Calculate activity rate for each node
}
\description{
Calculates the activity rate of the nodes (i.e., the number of times a node is active, i.e., ON, divided by the number of repeats) for a specified number of time-steps.
}
\usage{
calc_node_activities(net, method = c("BNp", "SDDS", "PEW"), params, steps,
                     repeats = 1000, initial_prob = NULL, last_step = FALSE,
                     asynchronous = TRUE, update_prob = NULL)
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
\item{steps}{The number of time-steps (non-negative integer) to simulate the network.}

\item{repeats}{The number of repeats (positive integer).}

\item{initial_prob}{The probability that each of the nodes is ON (1) in the initial state (time-step 0). If \code{NULL} (default), 0.5 is used as the probability for all nodes, meaning the initial state is randomly chosen based on a uniform distribution.}

\item{last_step}{If \code{TRUE}, only the node activity rates for the last time-step are returned. Otherwise, the entire activity rate trajectory of the nodes is returned.}

\item{asynchronous}{If \code{TRUE}, the asynchronous update scheme is used, where a single node is updated at each time-step. In this case, \code{update_prob} indicates update probabilities. If \code{FALSE}, the synchronous update scheme is utilized.}

\item{update_prob}{The probability of updating each variable (node) in each time-step when \code{asynchronous = TRUE}. If \code{NULL} (default), nodes are updated randomly based on a uniform distribution. If \code{asynchronous = FALSE}, this argument is ignored.}
}

\details{
By incorporating stochasticity into the update rule of a Boolean network and repeating the simulation several times, the average value of each node across repeats can be considered as a continuous variable. This approach transforms discrete binary variables into continuous ones, enabling continuous analysis methods applicable for studying the dynamic behavior of the Boolean network. This function calculates the average value (i.e., node activity) of each network node at each time-step.}

\value{
If \code{last_step = TRUE}, a vector with a length equal to the number of network nodes, representing the activity rate of each node at the last time-step, is returned. If \code{last_step = FALSE}, a matrix with \code{steps + 1} rows (where the first row corresponds to time-step 0) and \code{length(net$genes)} columns (representing node activities at each time-step) is returned. The order of the nodes in the vector or columns (depending on \code{last_step}) is the same as \code{net$genes}.
}

\references{
Golinelli, O., & Derrida, B. (1989). Barrier heights in the Kauffman model. \emph{Journal De Physique}, 50(13), 1587-1601.
Shmulevich, I., Dougherty, E. R., & Zhang, W. (2002). Gene perturbation and intervention in probabilistic Boolean networks. \emph{Bioinformatics}, 18(10), 1319-1331.

Shmulevich, I., Dougherty, E. R., & Zhang, W. (2002). Gene perturbation and intervention in probabilistic Boolean networks. \emph{Bioinformatics}, 18(10), 1319-1331.

Trairatphisan, P., Mizera, A., Pang, J., Tantar, A. A., Schneider, J., & Sauter, T. (2013). Recent development and biomedical applications of probabilistic Boolean networks. \emph{Cell communication and signaling}, 11, 1-25.

Murrugarra, D., Veliz-Cuba, A., Aguilar, B., Arat, S., & Laubenbacher, R. (2012). Modeling stochasticity and variability in gene regulatory networks. \emph{EURASIP Journal on Bioinformatics and Systems Biology}, 2012, 1-11.

Deritei, D., Kunšič, N., & Csermely, P. (2022). Probabilistic edge weights fine-tune Boolean network dynamics. \emph{PLoS Computational Biology}, 18(10), e1010536.
}

\examples{

# >>>>>>>>>>>>>>>>  Load network and generate random initial states  <<<<<<<<<<<<<<<<<

# Load the example network
data(lac_operon_net)

# Define plot function
plot_node_activities <- function(node_activities) {
  old_par <- par(no.readonly = TRUE)
  layout(matrix(c(1, 2), nrow = 1), width = c(4, 1))
  par(mar = c(5, 4, 4, 0))
  matplot(1:nrow(node_activities), node_activities, type = "l", frame = TRUE, lwd = 2,
      lty = 1, xlab = "Time-step", ylab = "Node activity")
  par(mar = c(5, 0, 4, 2))
  plot(c(0, 1), type = "n", axes = FALSE, xlab = "")
  legend("center", colnames(node_activities), col = seq_len(ncol(node_activities)),
      cex = 0.5, fill = seq_len(ncol(node_activities)))
  layout(1)
  par(old_par)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Method: BNp  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Define the parameters for the BNp method
params <- rep(0.05, length(lac_operon_net$genes))

# Get node activities after simulation using the BNp method
node_act <- calc_node_activities(lac_operon_net, method = "BNp", params = params,
    steps = 100, repeats = 10000)

# Plot node activities
plot_node_activities(node_act)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Method: SDDS  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Define the parameters for the SDDS method
props <- rep(0.95, length(lac_operon_net$genes))
params <- list(p00 = props, p01 = props, p10 = props, p11 = props)

# Get node activities after simulation using the SDDS method
node_act <- calc_node_activities(lac_operon_net, method = "SDDS", params = params,
    steps = 100, repeats = 10000)

# Plot node activities
plot_node_activities(node_act)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Method: PEW  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Extract edges from the network
edges <- extract_edges(lac_operon_net)

# Define the parameters for the PEW method
p_on <- runif(nrow(edges))
p_off <- runif(nrow(edges))
params <- list(p_on = p_on, p_off = p_off)

# Get node activities after simulation using the PEW method
node_act <- calc_node_activities(lac_operon_net, method = "PEW", params = params,
    steps = 100, repeats = 10000)

# Plot node activities
plot_node_activities(node_act)
}

