\name{pastboon-package}
\alias{pastboon-package}
\alias{pastboon}
\docType{package}
\title{
Simulation of Parameterized Stochastic Boolean Networks
}
\description{
Applying stochastic noise to Boolean networks is a useful approach for representing the effects of various perturbing stimuli on complex systems. A number of methods have been developed to control noise effects on Boolean networks using parameters integrated into the update rules. This package provides functions to simulate and analyze three such methods: \emph{Boolean network with perturbations} (BNp), described by Trairatphisan et al., \emph{stochastic discrete dynamical systems} (SDDS), proposed by Murrugarra et al., and \emph{Boolean network with probabilistic edge weights} (PEW), presented by Deritei et al. The package includes source code derived from the \code{BoolNet} package, which is licensed under the Artistic License 2.0.
}
\details{
Applying perturbations to a standard deterministic Boolean network involves altering its update rules. Manipulating the logical functions usually requires a thorough understanding of the reasoning behind the Boolean equations and may lead to a loss of the network's main functional characteristics, which often need to be preserved. An alternative approach to perturbing a Boolean network is to introduce stochastic noise and control its effect through a set of parameters integrated into the logical functions. This approach offers the advantage of allowing partial activation or inhibition of nodes.

In \code{pastboon}, three parameterization methods are implemented to control the stochastic noise effect on Boolean networks:

\itemize{
\item{BNp, \emph{Boolean network with perturbations} (Trairatphisan et al.)}
\item{SDDS, \emph{Stochastic discrete dynamical systems} (Murrugarra et al.)}
\item{PEW, \emph{Boolean network with probabilistic edge weights} (Deritei et al.)}
}

Given a Boolean network, its parameterization method, and the parameter values, useful insights can be gained from network simulations using the functions provided in this package. Node activity trajectories (the probability of each node being active at each time-step) can be calculated using \code{\link{calc_node_activities}}. Once node activity trajectories are available, the time-step at which the network reaches a steady-state distribution can be estimated using \code{\link{calc_convergence_time}}. Additionally, the states reached after starting a Boolean network from a given set of initial states can be sampled over specified time-steps using \code{\link{get_reached_states}}. The number of pairwise transitions between a given set of states can be obtained using \code{\link{count_pairwise_trans}}. Finally, the edges of a Boolean network can be extracted using \code{\link{extract_edges}}.

This package includes source code derived from the \code{BoolNet} package, which is licensed under the Artistic License 2.0. Specifically, the C code for simulating Boolean networks and its R interface code were initially taken from the \code{BoolNet} package but have been substantially altered (particularly the C code) to meet our purposes.
}

\author{
Mohammad Taheri-Ledari [aut, cre, cph] <mo.taheri@ut.ac.ir>

Kaveh Kavousi [ctb]

Sayed-Amir Marashi [ctb]

Authors of BoolNet [ctb]

Troy D. Hanson [ctb]
}

\references{
Trairatphisan, P., Mizera, A., Pang, J., Tantar, A. A., Schneider, J., & Sauter, T. (2013). Recent development and biomedical applications of probabilistic Boolean networks. \emph{Cell communication and signaling}, 11, 1-25.

Murrugarra, D., Veliz-Cuba, A., Aguilar, B., Arat, S., & Laubenbacher, R. (2012). Modeling stochasticity and variability in gene regulatory networks. \emph{EURASIP Journal on Bioinformatics and Systems Biology}, 2012, 1-11.

Deritei, D., Kunšič, N., & Csermely, P. (2022). Probabilistic edge weights fine-tune Boolean network dynamics. \emph{PLoS Computational Biology}, 18(10), e1010536.

Müssel, C., Hopfensitz, M., & Kestler, H. A. (2010). BoolNet—an R package for generation, reconstruction and analysis of Boolean networks. \emph{Bioinformatics}, 26(10), 1378-1380.
}
