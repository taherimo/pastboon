\name{lac_operon_net}
\alias{lac_operon_net}
\docType{data}
\title{
The lactose operon Boolean network
}
\description{
The *lactose* operon (*lac* operon) Boolean network as proposed by Veliz-Cuba and Stigler.
}
\usage{data(lac_operon_net)}
\details{
The data consists of a variable \code{lac_operon_net} of BoolNet's \code{BooleanNetwork} class describing *lac* operon gene regulary network with 10 genes and 3 inputs. The 3 inputs overally indicate the concentration of glucose and lactose. According to synchronous update, when extracellular glucose is available the *lac* operon is OFF (having one steady-state attractor where all genes are OFF), otherwise depending to extracellular lactose concentration the operon will be OFF, bistable (having two attractors), or ON (all genes are ON). The class \code{BooleanNetwork} is described in more detail in \code{\link{loadNetwork}}.
}
\source{
Veliz-Cuba, A., & Stigler, B. (2011). Boolean models can explain bistability in the lac operon. Journal of computational biology, 18(6), 783-794. [DOI:https://doi.org/10.1089/cmb.2011.0031](https://doi.org/10.1089/cmb.2011.0031)
}
\examples{
data(lac_operon_net)

# the network is stored in a variable called 'lac_operon_net'
print(lac_operon_net)

}
\keyword{E. coli
	lac operon
	Boolean network}
