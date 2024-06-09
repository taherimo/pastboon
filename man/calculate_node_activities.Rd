\name{calculate_node_activities}
\alias{calculate_node_activities}
\title{
Calculation of activity rate for each node
}
\description{
Calculates fraction of active states at each time step for each node among a number of repeats.
}
\usage{
calculate_node_activities(net, method=c("SDDS","BNp","PEW"),
                          params, steps, repeats=1000,
                          initial_prob= NULL, update_prob= NULL,
                          return_last_step= FALSE, asynchronous= TRUE)
}

\arguments{
\item{net}{A network structure of class \code{BooleanNetwork} from \pkg{BoolNet}.}

\item{method}{The parameterization method to be used. If "SDDS" is chosen, stochastic discrete dynamical system proposed by Murrugarra et al. will be used. If "BNp" is selected, Boolean network with perturbation initially presented by Golinelli and Derrida (and later by Shmulevich et al. and Trairatphisan et al.) will be utilized. If "PEW" is chosen, Boolean network with probabilistic edge weights proposed by Deritei et al. will be employed. For each of these cases a different pattern for the \code{params} argument is expected.}

\item{params}{The parameter values depending to the value of \code{method}. If \code{method= "SDDS"}, the \code{params} argument is expected to be a list of four equal-size vectors \code{p00}, \code{p01}, \code{p10}, and \code{p11}. If \code{method= "BNp"}, the \code{params} argument is expected to be vector of probabilities. Finally, if \code{params= "PEW"}, the \code{params} argument is expected to be a a list of two equall-size vectors \code{p_on} and \code{p_off}. The length of all probability vectors should be equal to the number of network nodes.}

\item{steps}{The number of time steps which should be a non-negatve integer.}

\item{repeats}{The number of time repeats which should be a positive integer.}

\item{initial_prob}{The probability that each of the variables are ON (1) in the initial state. If it is set to \code{NULL} (default value) then 0.5 is used as the probability for all nodes which means using uniform probability over state space of the network.}

\item{update_prob}{Probability of updating each variable in each time step is the case of \code{asynchronous= TRUE}. If it is set to \code{NULL} (default value) then nodes are randomly updated based on uniform distribution. If \code{asynchronous= FALSE} this argument is ignored.}

\item{return_last_step}{If \code{TRUE} only the node activity rates for the last time step is returned otherwise whole activity rate trajectory of the nodes are returned.}

\item{asynchronous}{If \code{TRUE}, the asynchronous update scheme is used where a single node is updated at each time-step. In this case \code{update_prob} indicates update probabuilities. If \code{FALSE}, synchronous update is utilized.}
}

\details{
In the current version of PARBONET, only \code{BooleanNetwork} networks are supported.
}

\value{
A matrix with \code{steps + 1} rows (the first row corresponds to initial node activities) and \code{length(net$genes)} columns that represents node activities in each time-step.
}

\references{
Murrugarra, D., Veliz-Cuba, A., Aguilar, B., Arat, S., & Laubenbacher, R. (2012). Modeling stochasticity and variability in gene regulatory networks. "EURASIP Journal on Bioinformatics and Systems Biology", 2012, 1-11.

Deritei, D., Kunšič, N., & Csermely, P. (2022). Probabilistic edge weights fine-tune Boolean network dynamics. "PLoS Computational Biology", 18(10), e1010536.

Golinelli, O., & Derrida, B. (1989). Barrier heights in the Kauffman model. "Journal De Physique", 50(13), 1587-1601.
Shmulevich, I., Dougherty, E. R., & Zhang, W. (2002). Gene perturbation and intervention in probabilistic Boolean networks. "Bioinformatics", 18(10), 1319-1331.

Trairatphisan, P., Mizera, A., Pang, J., Tantar, A. A., Schneider, J., & Sauter, T. (2013). Recent development and biomedical applications of probabilistic Boolean networks. "Cell communication and signaling", 11, 1-25.
}
\examples{

# load the example network
data(lac_operon_net)

###############################
# Example 1: The SDDS Method  #
###############################

# define the parameter for the SDDS method
props <- rep(0.95, length(lac_operon_net$genes)
params <- list(p00= props, p01= props, p10= props, p11= props)

# get node activities after simulation using the SDDS method and asynchronous update
node_act <- calculate_node_activities(lac_operon_net, method= "SDDS", params = params,
                            steps= 100, repeats = 10000)

# plot node activities
layout(matrix(c(1,2),nrow=1), width=c(4,1))
par(mar=c(5,4,4,0)) #No margin on the right side
matplot(1:nrow(node_act), node_act*100, type = "l", frame = TRUE, lwd = 2, lty = 1,
        xlab= "Time step", ylab= "Node activity (\%)")
par(mar=c(5,0,4,2)) #No margin on the left side
plot(c(0,1),type= "n", axes= FALSE, xlab= "Time step", ylab= "Node activity (\%)")
legend("center", colnames(node_act),col=seq_len(length(lac_operon_net$genes)),
        cex=0.5, fill=seq_len(length(lac_operon_net$genes)))

###############################
# Example 2: The BNp Method   #
###############################

# define the parameters for the BNp method
params <- rep(0.05, length(lac_operon_net$genes))

# get node activities after simulation using the BNp method and asynchronous update
node_act <- calculate_node_activities(lac_operon_net, method = "BNp", params = params,
                            steps=50, repeats = 5000)

# plot node activities

layout(matrix(c(1,2),nrow=1), width=c(4,1))

par(mar=c(5,4,4,0)) #No margin on the right side

matplot(1:nrow(node_act), node_act*100, type = "l", frame = TRUE, lwd = 2, lty = 1,
        xlab= "Time step", ylab= "Node activity (\%)")

par(mar=c(5,0,4,2)) #No margin on the left side

plot(c(0,1),type= "n", axes= FALSE, xlab= "Time step", ylab= "Node activity (\%)")

legend("center", colnames(node_act),col=seq_len(length(lac_operon_net$genes)),
    cex=0.5,fill=seq_len(length(lac_operon_net$genes)))


###############################
# Example 3: The PEW Method   #
###############################

# extract edges from the network
edges <- extract_edges(lac_operon_net)

# define the parameters for the PEW method
p_on <- runif(nrow(edges))
p_off <- runif(nrow(edges))
params <- list(p_on=p_on, p_off=p_off)

# get node activities after simulation using the PEW method and asynchronous update
node_act <- calculate_node_activities(lac_operon_net, method= "PEW", params = params,
                            asynchronous = TRUE, steps = 100, repeats = 10000)

# plot node activities

layout(matrix(c(1,2),nrow=1), width=c(4,1))

par(mar=c(5,4,4,0)) #No margin on the right side

matplot(1:nrow(node_act), node_act*100, type = "l", frame = TRUE, lwd = 2, lty = 1,
        xlab= "Time step", ylab= "Node activity (\%)")

par(mar=c(5,0,4,2)) #No margin on the left side

plot(c(0,1),type= "n", axes= FALSE, xlab= "Time step", ylab= "Node activity (\%)")

legend("center", colnames(node_act),col=seq_len(length(lac_operon_net$genes)),
       cex=0.5,fill=seq_len(length(lac_operon_net$genes)))

}
