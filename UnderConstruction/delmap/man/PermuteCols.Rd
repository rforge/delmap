\name{PermuteCols}
\alias{PermuteCols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Change Order of Marker Columns
}
\description{
\code{PermuteCols} column permutes the data in the object deldata.mat.
}
\usage{
PermuteCols(dmat= NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmat}{
An object of class \code{\link{delmap.data}}
}
}
\details{
\code{PermuteCols} permutes, randomly, the order of the columns of genotype data.
}
\value{
An object of class \code{\link{delmap.data}}.
}
\author{
Andrew W George
}

\examples{
# generate deletion data
Slist <-  SimDeletions( numlines = 50, plines = 0.2, nummarkers = 20,
    Enumdel = 4, p.missing = 0.02,  seed = 1)

# print the data, including the marker ordering
print(Slist[["nomissing"]])

# permute the data and print the data, including the new marker ordering
print(PermuteCols(Slist[["nomissing"]]))

}
