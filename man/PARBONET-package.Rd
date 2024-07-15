\name{PARBONET-package}
\alias{PARBONET-package}
\alias{PARBONET}
\docType{package}
\title{
Simulation of Parameterized Stochastic Boolean Networks
}
\description{
A useful approach to represent the effect of different perturbing stimuli on Boolean networks is stochastic noise. Some methods to apply stochastic noise on Boolean networks have been presented, which are based on controlling the noise effect using a set of parameters integrated into the update rules. In this package, some useful functions are provided to examine three such parameterization methods: \emph{Boolean network with perturbations} (BNp), described by Trairatphisan et al., \emph{stochastic discrete dynamical systems} (SDDS), proposed by Murrugarra et al., and \emph{Boolean network with probabilistic edge weights} (PEW), presented by Deritei et al. This package includes modified portions of source code from the \code{BoolNet} package.
}
\details{
\tabular{ll}{
Package: \tab PARBONET\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2024-10-01\cr
License: \tab Artistic-2.0\cr
}
Applying perturbations to a standard deterministic Boolean network is performed by altering its update rules. Manipulating the logical functions requires a proper understanding of the reasoning behind the Boolean equations. Moreover, it can lead to the loss of the main functional characteristics of the network, which usually need to be preserved. An alternative approach to perturb a Boolean network is to apply stochastic noise and control the noise effect through adjusting a set of parameters integrated into the logical functions. An important advantage of this approach is the possibility of implementing partial activation/inhibition of variables.

In \code{PARBONET}, three parameterization methods to control noise effect on Boolean networks are implemented:

\itemize{
\item{BNp, \emph{Boolean network with perturbations} (Trairatphisan et al.)}
\item{SDDS, \emph{Stochastic discrete dynamical systems} (Murrugarra et al.)}
\item{PEW, \emph{Boolean network with probabilistic edge weights} (Deritei et al.)}
}

Given a Boolean network, its parameterization method, and the values of the parameters, one can achieve useful information resulting from simulation of the network via the functions provided in this package. Node activity trajectories (the probability that a node is active at each time step) can be calculated using \code{\link{calc_node_activities}}. By having node activity trajectories, the time step at which a steady-state distribution is reached can be detected using \code{\link{calc_convergence_time}}. Also, the reached states after starting a Boolean network from a given set of initial states for some time steps can be sampled using \code{\link{get_reached_states}}. Moreover, the number of pairwise transitions between a given set of states can be obtained using \code{\link{count_pairwise_trans}}. Finally, the edges of a Boolean network can be extracted utilizing \code{\link{extract_edges}}.
}

\author{
Mohammad Taheri-Ledari [aut, cre, cph] <mo.taheri@ut.ac.ir>

Kaveh Kavousi [ctb]

Sayed-Amir Marashi [ctb]
}

\references{
Trairatphisan, P., Mizera, A., Pang, J., Tantar, A. A., Schneider, J., & Sauter, T. (2013). Recent development and biomedical applications of probabilistic Boolean networks. \emph{Cell communication and signaling}, 11, 1-25.

Murrugarra, D., Veliz-Cuba, A., Aguilar, B., Arat, S., & Laubenbacher, R. (2012). Modeling stochasticity and variability in gene regulatory networks. \emph{EURASIP Journal on Bioinformatics and Systems Biology}, 2012, 1-11.

Deritei, D., Kunšič, N., & Csermely, P. (2022). Probabilistic edge weights fine-tune Boolean network dynamics. \emph{PLoS Computational Biology}, 18(10), e1010536.

Müssel, C., Hopfensitz, M., & Kestler, H. A. (2010). BoolNet—an R package for generation, reconstruction and analysis of Boolean networks. \emph{Bioinformatics}, 26(10), 1378-1380.
}
