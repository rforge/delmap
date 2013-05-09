\name{dmappingchrm}
\alias{dmappingchrm}
\title{ Core routine dmapping}
\description{
\code{dmappingchrm} performs deletion mapping for a single chromosome. It is a core 
function in \code{\link{dmapping}}. 
}
\details{
Although this function provides low level ordering of marker loci on a single chromosome, 
\code{\link{dmapping}} is to be used for all deletion mapping (even when there may 
only be a single chromosome worth of data). 
}
\value{
A list object with elements \code{bestx}, \code{bestt},  
\code{besttvec},  and if reference ordering and block structure is provided, 
\code{bestdist}. 

Here, \itemize{
\item \code{bestx} is of object type \code{\link{dmdatachrm}}  and contains the best 
deletion map with imputed genotypes. 
\item \code{bestt} is a \code{numeric} value containing the (best) tour length for 
\code{bestx}.
\item \code{besttvec} is a \code{numeric} vector containing the best tour lengths at each 
iteration of the simulated annealing algorithm.
\item \code{bestdist} is a \code{\link{dmdist}} object for the (best) rearrangement distance 
for \code{bestx}.
}
}
\author{
Andrew W George
}
\seealso{
\code{\link{dmapping}}, \code{\link{dmdist}}
}
