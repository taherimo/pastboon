\name{PARBONET-package}
\alias{PARBONET-package}
\alias{PARBONET}
\docType{package}
\title{
Simulation of Parameterized Stochastic Boolean Networks
}
\description{
Provides functions for simulation of parameterized stochastic Boolean networks.
}
\details{
\tabular{ll}{
Package: \tab PARBONET\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2024-10-01\cr
License: \tab Artistic-2.0\cr
}
Adapting a Boolean network to different contexts finds many applications in different fields. Boolean networks usually are designed to describe average behaviour a system without taking into account its specificity in different conditions. For instance, the Boolean models which are presenred to portray molecular network of different cell types usually ignore heterogeneity within the cells. Altering such such Boolean networks to resemble specific cells is of great importance in precision medicine.

One way to manipulate dynamical behaviour of the Boolean networks is to incorporate stochasticity-related parameters in these networks and try to make the dynamical behaviour of the network as close as possible to our desired behavior by adjusting the parameter values. In the current version of \code{PARBONET}, three methods to parameterize the Boolean network are implemnted:

\itemize{
\item{SDDS (proposed by Murrugarra et al.)}
\item{BNp (proposed by Golinelli and Derrida)}
\item{PEW (proposed by Deritei et al.)}

}

Given a Boolean network, its parameterization method and the values of the parameters, one can achieve useful information resulting from simulation of the network via the functions provided in this package. Node activity trajectories (the probability that a node is active at each time-step) can be calculatd using \code{\link{calc_node_activities}}. The reached states after simulating a Boolean network for some pre-determined time steps can be sampled through \code{\link{get_reached_states}}. Pairwise state-transition frquencies for a set of states can be obtained by \code{\link{count_pairwise_transitions}}. By having node activity trajectories, the time-step that a steady-state distribution is reached can be detected using \code{\link{calc_convergence_time}}. Finally, the edges of a given Boolean network can be extracted by utilizing \code{\link{extract_edges}}.

}

\author{
Mohammad Taheri-Ledari

Maintainer: Mohammad Taheri-Ledari <mo.taheri@ut.ac.ir>
}

\references{
Murrugarra, D., Veliz-Cuba, A., Aguilar, B., Arat, S., & Laubenbacher, R. (2012). Modeling stochasticity and variability in gene regulatory networks. "EURASIP Journal on Bioinformatics and Systems Biology", 2012, 1-11.
\cr
\cr
Deritei, D., Kunšič, N., & Csermely, P. (2022). Probabilistic edge weights fine-tune Boolean network dynamics. "PLoS Computational Biology", 18(10), e1010536.
\cr
\cr
Golinelli, O., & Derrida, B. (1989). Barrier heights in the Kauffman model. "Journal De Physique", 50(13), 1587-1601.

Further references are listed in the corresponding reference sections.
}
