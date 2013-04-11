\name{plot.delmap.data}
\alias{plot.delmap.data}
\alias{plot}
\title{
Plot of Deletion Map
}
\description{
Produces a plot of the deletion map 
}
\usage{
\method{plot}{delmap.data}(x, main, \ldots)
}
\arguments{
  \item{x}{
an object of class 'delmap.dat' that has either been read in from file or simulated. 
}
  \item{main}{
title of plot.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
\code{plot.delmap.data} is a function for visualizing the deletion data as a deletion map. Data for the marker 
loci are plotted in column order. A solid dot denotes a deletion, a hollow dot denotes a missing value, and an empty cell 
denotes observed dna. A solid vertical line separates the marker blocks. 
}
\author{
James Totterdell
}
\seealso{
\code{\link{IdentifyMarkerBlocks}}
}
\examples{
# generate deletion data
Slist <-  SimDeletions( numlines = 50, plines = 0.2, nummarkers = 20,
    Enumdel = 4, p.missing = 0.02,  seed = 1)


## Clean data i.e. remove rows/columns without any deletions
Slist[["nomissing"]] <- CleanData(Slist[["nomissing"]])


# plot contents of data as a deletion map
plot(Slist[["nomissing"]])
}

