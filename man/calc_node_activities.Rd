\name{calc_node_activities}
\alias{calc_node_activities}
\title{
Calculate activity rate for each node
}
\description{
Calculates activity rate of the nodes (i.e. the number of times a node is active, i.e. ON, divided by the number of repeats) at each time-step.
}
\usage{
calc_node_activities(net, method = c("SDDS","BNp","PEW"), params, steps,
                     repeats = 1000, initial_prob = NULL, last_step = FALSE,
                     asynchronous = TRUE, update_prob = NULL)
}

\arguments{
\item{net}{A network structure of the class \code{BooleanNetwork} from the \pkg{BoolNet} package.}

\item{method}{The parameterization method to be used. For "SDDS", "BNp", amd
"PEW" the methods \emph{stochastic discrete dynamical systems}, \emph{Boolean network with perturbations}, and \emph{Boolean network with probabilistic edge weights} will be employed, respectively. For each of these cases a different pattern for the \code{params} argument is expected. For more information, see the Details section.
}

\item{params}{The parameter values depending on \code{method}. If \code{method = "SDDS"}, it is expected to be a list of four equal-size vectors of probabilities \code{p00}, \code{p01}, \code{p10}, and \code{p11} where all the vector lengths should be equal to the number of network nodes. If \code{method = "BNp"}, it should be a single vector of probabilities with the length equal to the number of network nodes. Finally, if \code{params = "PEW"}, it is required be a list of two equall-size vectors of probabilities \code{p_on} and \code{p_off} with the length equal to the number of network edges based on the order of the edges obtained from \code{\link{extract_edges}}.}

\item{steps}{The number of time-steps, which should be a non-negatve integer.}

\item{repeats}{The number of repeats, which should be a positive integer.}

\item{initial_prob}{The probability that each of the variables are ON (1) in the initial state (time-step 0). If \code{NULL} (default), 0.5 is used as the probability for all nodes which means the initial state is randomly chosen based on uniform distribution.}

\item{last_step}{If \code{TRUE}, only the node activity rates for the last time-step is returned otherwise whole activity rate trajectory of the nodes are returned.}

\item{asynchronous}{If \code{TRUE}, the asynchronous update scheme is used where a single node is updated at each time-step. In this case \code{update_prob} indicates update probabuilities. If \code{FALSE}, synchronous update is utilized.}

\item{update_prob}{Probability of updating each variable (node) in each time-step in the case of \code{asynchronous = TRUE}. If \code{NULL} (default) then nodes are randomly updated based on uniform distribution. If \code{asynchronous = FALSE} this argument is ignored.}
}

\details{
By incorporating any kind of stochasticity in the update rule of a Boolean network and repeating its simulation several times, one can consider the average value of each variable (node) across the repeats as a continuous variable. This way, the discrete binary variables are turned into continuous variables, which makes continuous analysis methods applicable to study the dynamic behavior of the Boolean network. This function calculates the average value (i.e. node activity) of each network node in each time-step.

Three Boolean network parameterization methods are implemented in this function: \emph{stochastic discrete dynamical systems} (SDDS) proposed by Murrugarra et al.; \emph{Boolean network with perturbations} (BNp) initially suggested by Golinelli and Derrida, later presented by Shmulevich et al. and Trairatphisan et al.; and \emph{Boolean network with probabilistic edge weights} (PEW) proposed by Deritei et al. All these methods actually different ways of including stochasticity in the Boolean networks.
}

\value{
If \code{last_step = TRUE}, a vector with the length equal to the number of network nodes, representing activity rate of each node at the last time-step. If \code{last_step= FALSE}, a matrix with \code{steps + 1} rows (the first row corresponds to the time-step 0) and \code{length(net$genes)} columns that represents node activities in each time-step. The order of the nodes in the columns is same as \code{net$genes}.
}

\references{
Murrugarra, D., Veliz-Cuba, A., Aguilar, B., Arat, S., & Laubenbacher, R. (2012). Modeling stochasticity and variability in gene regulatory networks. "EURASIP Journal on Bioinformatics and Systems Biology", 2012, 1-11.

Golinelli, O., & Derrida, B. (1989). Barrier heights in the Kauffman model. "Journal De Physique", 50(13), 1587-1601.
Shmulevich, I., Dougherty, E. R., & Zhang, W. (2002). Gene perturbation and intervention in probabilistic Boolean networks. \emph{Bioinformatics}, 18(10), 1319-1331.

Shmulevich, I., Dougherty, E. R., & Zhang, W. (2002). Gene perturbation and intervention in probabilistic Boolean networks. \emph{Bioinformatics}, 18(10), 1319-1331.

Trairatphisan, P., Mizera, A., Pang, J., Tantar, A. A., Schneider, J., & Sauter, T. (2013). Recent development and biomedical applications of probabilistic Boolean networks. \emph{Cell communication and signaling}, 11, 1-25.

Deritei, D., Kunšič, N., & Csermely, P. (2022). Probabilistic edge weights fine-tune Boolean network dynamics. \emph{PLoS Computational Biology}, 18(10), e1010536.
}

\examples{

# >>>>>>>>>>>>>>>>  Load network and generate random initial states  <<<<<<<<<<<<<<<<<

# Load the example network
data(lac_operon_net)

# Define plot function
plot_node_activities <- function(node_activities) {
  layout(matrix(c(1,2),nrow = 1), width = c(4,1))
  par(mar=c(5,4,4,0))
  matplot(1:nrow(node_activities), node_act, type = "l", frame = TRUE, lwd = 2,
  lty = 1, xlab= "time-step", ylab= "Node activity")
  par(mar = c(5,0,4,2))
  plot(c(0,1),type = "n", axes = FALSE, xlab = "")
  legend("center", colnames(node_act),col = seq_len(ncol(node_activities)),
      cex = 0.5, fill = seq_len(ncol(node_activities)))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Method: SDDS  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Define the parameter for the SDDS method
props <- rep(0.95, length(lac_operon_net$genes))
params <- list(p00= props, p01= props, p10= props, p11= props)

# Get node activities after simulation using the SDDS method
node_act <- calc_node_activities(lac_operon_net, method= "SDDS", params = params,
    steps = 100, repeats = 10000)

# Plot node activities
plot_node_activities(node_act)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Method: BNp  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Define the parameters for the BNp method
params <- rep(0.05, length(lac_operon_net$genes))

# Get node activities after simulation using the BNp method
node_act <- calc_node_activities(lac_operon_net, method = "BNp", params = params,
    steps=100, repeats = 10000)

# Plot node activities
plot_node_activities(node_act)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Method: PEW  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Extract edges from the network
edges <- extract_edges(lac_operon_net)

# Define the parameters for the PEW method
p_on <- runif(nrow(edges))
p_off <- runif(nrow(edges))
params <- list(p_on=p_on, p_off=p_off)

# Get node activities after simulation using the PEW method
node_act <- calc_node_activities(lac_operon_net, method= "PEW", params = params,
    steps = 100, repeats = 10000)

# Plot node activities
plot_node_activities(node_act)
}
