\name{plot.dmdata}
\alias{plot.dmdata}
\alias{plot}
\title{
Colour Plot of Deletion Map
}
\description{
Produces a colour plot of the deletion map 
}
\usage{
\method{plot}{dmdata}(x, main, \ldots)
}
\arguments{
  \item{x}{
an object of class \code{dmdata} that is allowed to have missing genotypes.
}
  \item{main}{
title of plot.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
A deletion map plot is produced from \code{x}. If \code{x} contains missing genotypes (NA), 
then a simple black and white plot is produced where solid grey cirlces denote observed 
deletions and unfilled grey circles denote missing genotypes. 


If \code{x} does not contain any missing genotypes, then a colour plot is produced. Markers that do not vary in their genotype across the lines
have marker labels of the same colour. A solid grey circle denotes an observed deletion. A solid
red circle denotes an imputed deletion. A unfilled grey circle denotes an unobserved genotype that
has not been imputed. A white square denotes an observed non-deletion. A solid red square
denotes an imputed non-deletion.  A solid vertical line separates the marker blocks.


Both the black and white plot and colour plot assumed the marker order to be the column order of 
\code{x}. 


}
\author{
Andrew George
}
\seealso{
\code{\link{IdentifyMarkerBlocks}}
}
\examples{
# generate map for a single chromosome that has 20 equally distance marker loci
map <- sim.dmmap(100, 20)

# generate deletion data
Slist <-  sim.dmdata( numlines = 50, plines = 0.2, map = map, 
    Enumdel = 4, p.missing = 0.02,  seed = 1)


## Clean data i.e. remove rows/columns without any deletions
Slist[["nomissing"]] <- CleanData(Slist[["nomissing"]])


# this produces a colour plot of the data
plot(Slist[["nomissing"]])

# this produces a black & white plot of the data
plot(Slist[["missing"]])

}

