\name{sim.dmdata}
\alias{sim.dmdata}
\title{
Simulation of deletion data
}
\description{
\code{sim.dmdata} generates deletion data according to the linkage map \code{map}. If \code{map} 
is integer, then all marker loci are assumed to be on the same chromosome. Data are 
generated where the rows correspond to the plants and the columns correspond to the marker loci.
}
\usage{
sim.dmdata(numlines = NULL, plines = 0.5, map = NULL, numblocks=1, 
           Enumdel = 4, p.missing = 0, seed = NULL)
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
a list object that specifies the marker map. The marker map may be read in from 
file with \code{\link{read.files}} or simulated with 
\code{\link{sim.dmmap}}. If an integer number, the marker loci are  assumed to
belong to the same chromosome. 
}
  \item{numblocks}{
the minimum number of blocks (See \code{\link{idmarkerblocks}}) belonging to
a chromosome. }
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
Data are generated either under an assumed marker map, a chromosome at a time or where the marker loci 
are assumed to belong to the same chromosome.  The marker 
map is either read from file with \code{\link{read.files}} or simulated with 
\code{\link{sim.dmmap}}. If the map is of no importance, then assign \code{map} an 
integer number which corresponds to the  number of marker loci for which data are to be simulated. 

\code{numblocks} is the minumum number of marker blocks per chromosome. 
Marker blocks occur when there are no overlapping deletions across the lines. 
If the probability of a deletion is low and/or the length of a deletion is small, then a 
chromosome many contain many more blocks than \code{numblocks}. 

The length of a deletion is random and follows a 
poisson distribution with mean \code{Enumdel}. If \code{p.missing} is set, then the 
data will contain missing genotypes (\code{NA}). 
}
\value{
a list is returned that contains a single element \code{nomissing} if  \code{p.missing} is 0
or two elements \code{nomissing} and \code{missing} if \code{p.missing} is non-zero. 

Each element is an \code{\link{dmdata}} object which is itself a \code{list} object with 
the number of elements corresponding to the number of chromosomes in the marker map \code{map} and 
an element for the marker map. 

}
\author{
Andrew W. George
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{sim.dmmap}}
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

# generate deletion data
# a list is created with elements "missing" and "nomissing"
#corresponding to deletion data with and without missing genotypes,
# respectively.
S <- sim.dmdata( numlines = 500, plines = 0.2, map = map, 
    numblocks=2, Enumdel = 8, p.missing = 0.02,  seed = 11)

 numlines = 500; plines = 0.2; map = map
    numblocks=2; Enumdel = 8; p.missing = 0.02;  seed = 11


# print contents of matrix
print(S$nomissing)

# print contents of matrix
print(S$missing)



}
 
