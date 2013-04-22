\name{print.dmdata}
\alias{print.dmdata}
\alias{print}
\title{
Print Deletion Data
}
\description{
'print' prints the deletion data.
}
\usage{
\method{print}{dmdata}(x, \ldots)
}
\arguments{
  \item{x}{
an object of class 'delmap.dat' that has either been read in from file or simulated. 
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
'print.dmdata' is a customized print function for printing the contents of the 'dmdata' object. The first
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
S <- sim.dmdata( numlines = 500, plines = 0.2, map = map, 
    Enumdel = 8, p.missing = 0.02,  seed = NULL)

# print contents of data where all genotypes are observed
print(S[["nomissing"]])

#print contents of data where there are missing genotypes
print(S[["missing"]])

}

