\name{print.dmdatachrm}
\alias{print.dmdatachrm}
\title{
Print Deletion Data
}
\description{
Prints the deletion data assuming all marker loci belong to the same 
chromosome.
}
\usage{
\method{print}{dmdatachrm}(x, \ldots)
}
\arguments{
  \item{x}{
an object of class \code{\link{dmdatachrm}}. 
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
\code{print.dmdatachrm} is a customized print function for printing 
the contents of the \code{\link{dmdatachrm}} object. The first
 five lines of data, plant names, and marker names are displayed for 
each chromosome. 
}
\author{
Andrew W. George
}
\seealso{
\code{\link{summary}}
}
\examples{
# generate map for a single chromosome with 20 equally distant marker loci
map <- sim.dmmap(rep(100,5),20)

# generate deletion data
# a list is created with elements "missing" and "nomissing"
#corresponding to deletion data with and without missing genotypes,
# respectively.
S <- sim.dmdata( numlines = 500, plines = 0.2, map = map, 
    Enumdel = 8, p.missing = 0.02,  seed = NULL)

# print the contents of chromosome 5 say that have missing genotypes
class(S[["missing"]][[5]])
print(S[["missing"]][[5]])


# to print the contents of all the chromosomess
print(S[["missing"]])



}

