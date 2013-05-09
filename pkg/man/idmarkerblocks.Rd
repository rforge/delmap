\name{idmarkerblocks}
\alias{idmarkerblocks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify Marker Blocks
}
\description{
Identifies all marker blocks in  the deletion data. The deletion 
data cannot contain any missing genotypes.
}
\usage{
idmarkerblocks(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{\link{dmdata}} or \code{\link{dmdatachrm}}  which contains the genotype data 
in 0/1 form and has no missing genotypes.
}
}
\details{
A marker block is a group of marker loci that can be ordered relative to each other. 
Marker blocks occur when there are no overlapping deletions that span the marker loci in different 
blocks. Marker blocks and be rearranged in any order without affecting the goodness of the ordering.   


There may be multiple marker blocks for a set of deletion data. The markers within the blocks can be ordered but the marker loci across blocks cannot be ordered relative to each other. 

Missing genotypes are not allowed and must first be imputed with \code{\link{impute}}
}
\value{
If \code{x} is of class \code{\link{dmdatachrm}}, then a named integer vector of length the number of marker loci 
is returned. 

If \code{x} is of class \code{\link{dmdata}}, then a list object is returned where the number of elements in the
list is the number of chromosomes for which there are deletion data. Each element in the list is a named 
integer vector of length the number of marker loci on the chromosome. 

As an example, a vector of 1,1,1,2,2,3,3 means that markers 1,2, and 3 belong to a block, markers 
4 and 5 belong to a seperate marker block, and markers 6 and 7 belong to another block of markers.  
}

\author{
Andrew W. George
}
\note{
The marker blocks may change if marker genotypes are being imputed.
}
\seealso{
\code{\link{idmarkersuborder}}
}
\examples{
# generate map with 10 chromosomes, each of length 
# 100 cM and each with 50  equally spaced marker loci. 
map <- sim.dmmap(rep(100,10), 50)
# generate deletion data
Slist <-  sim.dmdata( numlines = 50, plines = 0.2, map = map,
    Enumdel = 4, p.missing = 0.02,  seed = 1)

# class of object that does not contain any missing data
class(Slist[["nomissing"]])

## Clean data i.e. remove rows/columns without deletion data
Slist[["nomissing"]] <- cleandata(Slist[["nomissing"]])

# print the block structure for each chromosome
idmarkerblocks(Slist[["nomissing"]])
}
