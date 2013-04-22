\name{sim.dmdata}
\alias{sim.dmdata}
\title{
Simulation of deletion data
}
\description{
\code{sim.dmdata} returns a matrix of genotype data where the rows correspond to 
plants and the columns correspond to arker loci. A 
column (i.e., marker locus) can contain genotypes 0/1/NA where 0 is no 
deletion, 1 is a deletion, and NA is missing.   
}
\usage{
sim.dmdata(numlines = NULL, plines = 0.5, map = NULL, numblocks=1, Enumdel = 4, p.missing = NULL, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numlines}{
number of inbred lines/plants.
}
  \item{plines}{
the probability of a plant carrying a deletion. 
}
  \item{map}{
a list object that specifies the marker map. If map is unknown, use \code{\link{sim.dmmap}} to 
create a generic map.
}
  \item{numblocks}{
the minimum number of blocks (See \code{\link{IdentifyMarkerBlocks}}) in which each 
chromosome will have. }
  \item{Enumdel}{
the average length (in terms of number of marker loci) for a deletion.
}
  \item{p.missing}{
the probability of a genotype being unobserved.
}
  \item{seed}{
an integer number that if set, will produce reproducable replicates. 
}
}
\details{
Data are generated under an assumed marker map, a chromosome at a time.  The marker 
map is either read from file with \code{\link{read.dmmap}} or simulated with 
\code{\link{sim.dmmap}}. If the map is of no importance, then 
\code{sim.dmmap(100,x)} can be used to create a simple single chromosome map with   \code{x} 
marker loci.  

\code{numblocks} is the minumum number of marker blocks per chromosome. 
Marker blocks occur when there are no overlapping deletions across the lines. 
If the probability of a deletion is low and/or the length of a deletion is small, then a 
chromosome many contain many more blocks than \code{numblocks}. 

The length of a deletion is random and follows a 
poisson distribution with mean \code{Enumdel}. If \code{p.missing} is set, then the 
data will contain missing genotypes (\code{NA}). 
}
\value{
a list is returned that contains only a single element if \code{p.missing} is \code{NULL} 
or two elements if \code{p.missing} is set. The elements are:
\item{nomissing}{is a matrix of genotype data 0/1 where 0 is no deletion and
1 is a deletion. The number of rows corresponds to the number of plants 
in the simulated population and the number of columns is the number of 
genetic marker loci. This will be the only element in output list if 
p.missing is \code{NULL}.}
\item{missing}{is a matrix of genotype data 0/1 that  
is equivalent to "nomissing" but where marker genotypes have 
been replaced, randomly, and set to NA.} 
}
\author{
Andrew W. George
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{plot}}
}
\examples{

# To generate a inbred population of 500 plants 
# where genotype data are collected on 10 chromosomes 
# where there are 100 markers per chromosome. Each chromosome
# has length 250 cM. 
# The probability of a plant having a genomic deletion is 0.2, 
# the average length of a deletion is 8 marker loci, and 
# the probability of a marker genotype being missing is 0.02. Each
# chromosome contains at least 2 blocks. 
# use:
map <- sim.dmmap(rep(250,10),100)

S <- sim.dmdata( numlines = 500, plines = 0.2, map = map, 
    numblocks=2, Enumdel = 8, p.missing = 0.02,  seed = 11)

# print contents of matrix
print(S$nomissing)

# print contents of matrix
print(S$missing)

# plot contents
plot(S$nomissing)



}
 
