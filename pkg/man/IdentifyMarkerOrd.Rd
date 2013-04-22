\name{IdentifyMarkerOrd}
\alias{IdentifyMarkerOrd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify Marker Orderings
}
\description{
This function identifies the different possible marker orders that are possible with a set of deletion data and an assumed 
marker order.  The assumed marker order is derived from the order of the marker columns in the object \code{dmat}. 
}
\usage{
IdentifyMarkerOrd(dmat = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmat}{
An object of class \code{\link{dmdata}} which contains the genotype data in 0/1 form and has no missing genotypes.
}
}
\details{
In some cases, it is possible to change the ordering of the marker loci (ie the order of the columns of genotype data) without changing the overall pattern of deletions.  With \code{IdentifyMarkerOrd} we are able to find all possible orderings that are consistent with an assumed marker ordering. 

The different marker orderings may occur in marker blocks. These blocks are identified with the function \code{IdentifyMarkerBlocks}.  Marker blocks are groupings of marker loci that, according to the deletion data, can occur in any order. \code{IdentifyMarkerOrd} orders marker loci within marker blocks.

Missing genotypes are not allowed and must first be imputed with \code{\link{ImputeMissingGeno} }
}
\value{
A list is returned with entries
\item{blocks}{is an integer vector which gives which block each marker loci belongs}
\item{orders}{is a list where the number of entries corresponds to the number of marker blocks. Each list entry 
is an integer vector where the number of integer elements corresponds to the number of marker loci in the marker block. 
The integer elements indicate which marker loci can be permuted without chaning the overall pattern of deletions (given the assumed 
marker ordering). For example, a vector of 1,1,2,3,3,3 means that markers 1 and 2 can be swapped without affecting the pattern of 
deletions, marker 3 cannot be swapped with any other marker, and markers 4, 5, and 6 can be permuted.   
}
}
\author{
Andrew W. George
}
\note{
The marker blocks may change if marker genotypes are being imputed.
}
\seealso{
\code{\link{ImputeMissingGeno}}; \code{\link{IdentifyMarkerBlocks}}
}
\examples{
# generate map for a single chromosome with 20 equally distant marker loci
map <- sim.dmmap(100,20)

# generate deletion data
Slist <-  sim.dmdata( numlines = 50, plines = 0.2, map = map,
    Enumdel = 4, p.missing = 0.02,  seed = 1)

# class of object that does not contain any missing data
class(Slist[["nomissing"]])

## Clean data i.e. remove rows/columns without deletion data
Slist[["nomissing"]] <- CleanData(Slist[["nomissing"]])

# print the block and possible marker orderings with each block for
# the simulated deletion data
IdentifyMarkerOrd(Slist[["nomissing"]])
}
