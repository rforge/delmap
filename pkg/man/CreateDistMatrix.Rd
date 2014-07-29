\name{CreateDistMatrix}
\alias{CreateDistMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distance Measure Computation
}
\description{
This function computes the Manhattan distance matrix from 0/1 deletion data 
for a given chromosome. It is assumed the
deletion data does not contain any missing genotypes. 
}
\usage{
CreateDistMatrix(x = NULL, dmethod="manhattan")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a \code{\link{dmdatachrm}} object or matrix that 
contains only 0's and 1's for a single chromosome.
}
   \item{dmethod}{
the distance measure to be used. This must be one of
          "euclidean", "maximum", "manhattan", "canberra",
          "binary" or "minkowski".  Any unambiguous substring can
          be given.
}
}
\details{
Calculates the Manhattan distance between all pairs of marker loci 
(columns) on a given chromosome. The function 
will produce an error if the data contains missing values. 
}
\value{
Returns an object of class 'dist'.
}
\author{
Andrew W. George
}

\examples{
## Example 1
## Simple example of where deletion data are contained in 
## a matrix object. 

# deletion data where 1 is a deletion
d <- matrix(data=c(0,0,1,1, 
                        0,0,0,1, 
                        1,0,0,0, 
                        0,1,0,0,  
                        1,1,1,0), nrow=5)
# convert  into dmdatachrm  object
dmp <- as.dmdatachrm(d)

# calculate distance matrix for deletion data
CreateDistMatrix(dmp)


## Example 2
## A more complicated example where deletion data are recorded for 
## 100 marker loci across 10 chromosome. 

# generate map
map <- sim.dmmap(rep(250,10), 10)

# generate deletion data
S <-  sim.dmdata( numlines = 500, plines = 0.2, map = map,
    numblocks=2, Enumdel = 8, p.missing = 0.02,  seed = 11)

# calculate distance for full data (S[["nomissing"]]) 
# for chromosome 2 say
CreateDistMatrix(S[["nomissing"]][[2]])



}

