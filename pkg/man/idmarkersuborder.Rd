\name{idmarkersuborder}
\alias{idmarkersuborder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify Marker Orderings
}
\description{
Identifies the different possible marker orders that are possible with a set of deletion data and an assumed 
marker order.  The assumed marker order is derived from the order of the 
marker columns in the object \code{x}.  
}
\usage{
idmarkersuborder(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{\link{dmdata}}  or \code{\link{dmdatachrm}}.
}
}
\details{
It will often be the case where a dervied ordering for the marker loci 
will not in fact be unique. This lack of uniqueness occurs when 
the marker loci are in blocks  (see \code{\link{id}})  and when 
marker loci share identical marker data within the blocks. 

\code{idmarkersuborder} gives the marker ordering within each block across 
chromsomes. For example, supposed the best ordering of marker loci 
M1, \ldots, M7 is M1 M2 M3 M4 M5 M6 M7 where these loci belong to the same 
marker block. Lets assume that the marker data for M1 and M2 are the 
same. Also, lets assume the marker data for M5 and M6 are the same. 
Then the output from \code{idmarkersuborder} would be 1 1 2 3 4 4 5. 

Missing genotypes are not allowed and must first be imputed with \code{\link{ImputeMissingGeno} }
}
\value{
If \code{x} is of class \code{\link{dmdatachrm}} then
a list is returned with entries:
\itemize{
\item{blocks}{is an integer vector which gives which block each marker loci belongs}
\item{orders}{is a list where the number of entries corresponds to the number of marker blocks. Each list entry 
is an integer vector where the number of integer elements corresponds to the number of marker loci in the marker block. 
The integer elements indicate which marker loci can be permuted without chaning the overall pattern of deletions (given the assumed 
marker ordering). For example, a vector of 1,1,2,3,3,3 means that markers 1 and 2 can be swapped without affecting the pattern of 
deletions, marker 3 cannot be swapped with any other marker, and markers 4, 5, and 6 can be permuted.   
}}

If \code{x} is of class \code{\link{dmdata}}, then a list of a list 
is returned. Each element in the list corresponds to a chromosome and 
the elements are themselves a list as described above. 
}
\author{
Andrew W. George
}
\note{
The marker blocks and marker orderings 
may change if marker genotypes are being imputed.
}
\seealso{
\code{\link{ImputeMissingGeno}}; \code{\link{idmarkerblocks}}
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

# print the marker ordering within each block for each chromosome
idmarkersuborder(Slist[["nomissing"]])
}
