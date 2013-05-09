\name{print.dmdata}
\alias{print.dmdata}
\alias{print}
\title{
Print Deletion Data
}
\description{
Prints the deletion data for each chromosome.
}
\usage{
\method{print}{dmdata}(x, \ldots)
}
\arguments{
  \item{x}{
an object of class \code{\link{dmdata}}.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
\code{print.dmdata} is a customized print function for printing 
the contents of the \code{\link{dmdata}} object. For each 
chromosome, 
the first
 five lines of data, plant names, and marker names are displayed. 
}
\author{
Andrew W. George
}
\seealso{
\code{\link{summary}}
}
\examples{
# generate map for a single chromosome with 20 equally distant marker loci
map <- sim.dmmap(100,20)

# generate deletion data
# a list is created with elements "missing" and "nomissing"
#corresponding to deletion data with and without missing genotypes,
# respectively.
S <- sim.dmdata( numlines = 500, plines = 0.2, map = map, 
    Enumdel = 8, p.missing = 0.02,  seed = NULL)

# print contents of data where all genotypes are observed
print(S[["nomissing"]])

#print contents of data where there are missing genotypes
print(S[["missing"]])

}

