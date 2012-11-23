\name{plot.delmap.data}
\alias{plot.delmap.data}
\alias{plot}
\title{
Plot of Deletion Map
}
\description{
Produces a coloured plot of the deletion map 
}
\usage{
\method{plot}{delmap.data}(x, \ldots)
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
\code{plot.delmap.data} is a function for visualizing the deletion data as a deletion map. Data for the marker 
loci are plotted in column order. A different colour scheme is used for each separate marker block. 
}
\author{
Andrew W. George
}
\seealso{
\code{\link{IdentifyMarkerBlocks}}
}
\examples{
# generate deletion data
S <- SimDeletions( numlines = 500, plines = 0.2, nummarkers = 20,
    Enumdel = 8, p.missing = 0.02,  seed = NULL)

# plot contents of data as a deletion map
plot(S$nomissing)
}

