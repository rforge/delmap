\name{subsetdata}
\alias{subsetdata}
\title{
Subset Deletion Data
}
\description{
Returns a dmdata object subsetted on chromosomes, marker loci, and/or rows.
}
\usage{
subsetdata(x=NULL, keep.chrm=NULL, drop.mrks=NULL, drop.rows=NULL)
}
\arguments{
  \item{x}{
an object of class \code{\link{dmdata}}. 
}
 \item{keep.chrm}{
numeric, integer, or character values that references the chromosomes to be kept.
}
 \item{drop.mrks}{
numeric, integer, or character values that reference the marker loci to be removed.
}
 \item{drop.rows}{
numeric, integer, or character values that reference the rows to be removed.
}
}
\details{
The function is used to form a subset of an object of class \code{\link{dmdata}}.
Data can be subsetted on chromosomes. Marker and/or rows can be removed. 

If marker (columns) and/or rows are being removed, \code{keep.chrm} must be specified. Here,  only 
a single chromosome is allowed to be specified. 

}
\value{
An object of class \code{\link{dmdata}} but for the specified chromosome(s) and with possibly rows and columns removed. The map 
is adjusted accordingly, along with the \code{nchrm} attribute. 
}
\author{
Andrew W. George
}
\seealso{
\code{\link{cleandata}}
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

# form a new dmdata object for the fifth chromosome
Ssub <- subsetdata(S[["missing"]], keep.chrm=5)
Ssub


# now remove some marker loci and some rows of data and adjust the map
Ssub <- subsetdata(Ssub, keep.chrm=1, drop.mrks=1:10, drop.rows=c(2,6,29))
Ssub


}

