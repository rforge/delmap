\name{sim.map}
\alias{sim.map}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate marker map
}
\description{
Simulate the positions of marker loci on a marker map
}
\usage{
sim.map(len=rep(100,20), n.mar=10, eq.spacing=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{len}{
a vector specifying the chromosome lengths (in cM).
}
  \item{n.mar}{
a vector specifying the number of markers per chromosome.
}
  \item{eq.spacing}{
a vector specifying the number of markers per chromosome.
}
}
\details{
Marker positions are simulated as iid Uniform(0,L). 
If \code{len} or \code{n.mar} has just one element, it is expanded to the length of the other argument. 
If they both have just one element, only one chromosome is simulated.

If \code{eq.spacing} is \code{TRUE}, markers are equally spaced between 0 and L. 


}
\author{
modified by Andrew George. Orginal version by Karl W. Broman
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{sim.dmdata}}
}
\examples{
# simulate 10 chromosomes of differing lengths but each with 10 
# marker loci
map <- sim.map(c(100,90,80,40,40,50,120,200,150,180), 10)

map



}
 
