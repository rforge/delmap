\name{print.delmap.data}
\alias{print.delmap.data}
\alias{print}
\title{
Print Deletion Data
}
\description{
'print' prints the deletion data.
}
\usage{
\method{print}{delmap.data}(x, \ldots)
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
'print.delmap.data' is a customized print function for printing the contents of the 'delmap.data' object. The first
 five lines of data, plant names, and marker names are displayed. 
}
\author{
Andrew W. George
}
\seealso{
\code{\link{summary}}
}
\examples{
# generate deletion data
S <- SimDeletions( numlines = 500, plines = 0.2, nummarkers = 20,
    Enumdel = 8, p.missing = 0.02,  seed = NULL)

# print contents of data where all genotypes are observed
print(S[["nomissing"]])

#print contents of data where there are missing genotypes
print(S[["missing"]])

}

