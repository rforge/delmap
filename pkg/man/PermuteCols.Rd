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
An object of class \code{\link{dmdata}}
}
}
\details{
\code{PermuteCols} permutes, randomly, the order of the columns of genotype data.
}
\value{
An object of class \code{\link{dmdata}}.
}
\author{
Andrew W George
}

\examples{
# generate map for a single chromosome with 20 equally spaced marker loci
map <- sim.dmmap(100, 20)

# generate deletion data
Slist <-  sim.dmdata( numlines = 50, plines = 0.2, map = map, 
    Enumdel = 4, p.missing = 0.02,  seed = 1)

# print the data, including the marker ordering
print(Slist[["nomissing"]])

# permute the data and print the data, including the new marker ordering
print(PermuteCols(Slist[["nomissing"]]))

}
