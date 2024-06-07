\name{PARBONET-package}
\alias{PARBONET-package}
\alias{PARBONET}
\docType{package}
\title{
Simulation of Parameterized Stochastic Boolean Networks
}
\description{
Functions for simulation of parameterized stochastic Boolean networks.
}
\details{
\tabular{ll}{
Package: \tab PARBONET\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2024-10-01\cr
License: \tab Artistic-2.0\cr
}
Adapting a Boolean networks to different contexts finds many applications in different fields. Boolean networks usually describe a general phenomena regardless of its specificity in different conditions. For instance Boolean models of molecular networks usually describe behavour a known cell type by ignoring its heterogeneity among different cases and also its heterogeneity in a tissue. Altering such a Boolean network to adapt it to a specific condition (e.g. a patient or a cell sub-type) is of great importance in precision medicine.

One way to manipulate dynamical behaviour of the Boolean networks is to incorporate stochasticity-related parameters in the Boolean networks and trying to make the dynamical behaviour of the system as close as possible to our desired behavior by setting the parameter values. In the current version of this package, three methods of Boolean network parameterization are implemnted which include:

\itemize{
\item{SDDS (proposed by Murrugarra et al.)}
\item{BNp (proposed by Golinelli and Derrida)}
\item{PEW (proposed by Deritei et al.)}

}

One can observe node activity rates for an ensemble of the Boolean networks with the same parameter values using the \cite{\link{get_node_activities}} function. Additionally, one can sample reached states after simulating a Boolean network for some pre-determined time steps using the \cite{\link{get_reached_states}} function.

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
