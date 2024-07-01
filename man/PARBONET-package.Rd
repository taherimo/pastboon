\name{PARBONET-package}
\alias{PARBONET-package}
\alias{PARBONET}
\docType{package}
\title{
Simulation of Parameterized Stochastic Boolean Networks
}
\description{
Provides functions for simulation of parameterized stochastic Boolean networks. This package includes modified portions of source code from the \pkg{BoolNet} package.
}
\details{
\tabular{ll}{
Package: \tab PARBONET\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2024-10-01\cr
License: \tab Artistic-2.0\cr
}
Adapting a Boolean network to different contexts finds many applications in different fields. Boolean networks usually are designed to describe average behaviour a system without taking into account its specificity in different conditions. For instance, the Boolean models which are presenred to portray molecular network of different cell types usually ignore heterogeneity within the cells. Altering such such Boolean networks to resemble specific cell types is of great importance in precision medicine.

One way to manipulate dynamical behaviour of the Boolean networks is to incorporate stochasticity-related parameters in these networks and try to make the dynamical behaviour of the network as close as possible to our desired behavior by adjusting the parameter values. In the current version of \code{PARBONET}, three methods to parameterize the Boolean network are implemnted:

\itemize{
\item{BNp, \emph{Boolean network with perturbations} (Trairatphisan et al.)}
\item{SDDS, \emph{Stochastic discrete dynamical systems} (Murrugarra et al.)}
\item{PEW, \emph{Boolean network with probabilistic edge weights} (Deritei et al.)}

}

Given a Boolean network, its parameterization method and the values of the parameters, one can achieve useful information resulting from simulation of the network via the functions provided in this package. Node activity trajectories (the probability that a node is active at each time-step) can be calculatd using \code{\link{calc_node_activities}}. The reached states after simulating a Boolean network for some pre-determined time steps can be sampled through \code{\link{get_reached_states}}. Pairwise state-transition frquencies for a set of states can be obtained by \code{\link{count_pairwise_trans}}. By having node activity trajectories, the time-step that a steady-state distribution is reached can be detected using \code{\link{calc_convergence_time}}. Finally, the edges of a given Boolean network can be extracted by utilizing \code{\link{extract_edges}}.

}

\author{
Mohammad Taheri-Ledari

Maintainer: Mohammad Taheri-Ledari <mo.taheri@ut.ac.ir>
}

\references{
Müssel, C., Hopfensitz, M., & Kestler, H. A. (2010). BoolNet—an R package for generation, reconstruction and analysis of Boolean networks. \emph{Bioinformatics}, 26(10), 1378-1380.

Trairatphisan, P., Mizera, A., Pang, J., Tantar, A. A., Schneider, J., & Sauter, T. (2013). Recent development and biomedical applications of probabilistic Boolean networks. \emph{Cell communication and signaling}, 11, 1-25.

Murrugarra, D., Veliz-Cuba, A., Aguilar, B., Arat, S., & Laubenbacher, R. (2012). Modeling stochasticity and variability in gene regulatory networks. \emph{EURASIP Journal on Bioinformatics and Systems Biology}, 2012, 1-11.

Deritei, D., Kunšič, N., & Csermely, P. (2022). Probabilistic edge weights fine-tune Boolean network dynamics. \emph{PLoS Computational Biology}, 18(10), e1010536.
}
