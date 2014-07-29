\name{dmplot}
\alias{dmplot}
\title{
Colour Plot of the Deletion Map
}
\description{
Produces a colour plot of the deletion map. 
}
\usage{
dmplot(x=NULL, keep.chrm=NULL, drop.mrks=NULL, drop.rows=NULL, keep.mrks=NULL, keep.rows=NULL, sumx=NULL, main=NULL, cx=0.5, index=NULL, ...)
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

 \item{keep.mrks}{
numeric, integer, or character values that reference the marker loci to be kept.
}
 \item{keep.rows}{
numeric, integer, or character values that reference the rows to be kept.
}
  \item{sumx}{the X axis labels to be used. Must be one of \code{NULL}, \code{cM}, \code{NAs}, \code{dels}.
}
  \item{main}{
title of plot.
}
  \item{cx}{
size of dots. Larger the value, the larger the dots. 
}
  \item{index}{
if specified, a pdf file called Rplots\code{index}.pdf produced. This index is useful when a large number of plots are to be produced.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
NEEDS WORK

A deletion map plot is produced from \code{x}. If \code{x} contains missing genotypes (NA), 
then a  plot is produced where solid grey cirlces denote observed 
deletions and blue circles denote missing genotypes. 

If \code{x} does not contain any missing genotypes, then the marker labels are colour coded. Markers 
that do not vary in their genotype across the lines
have marker labels of the same colour.  Here,  grey circles denotes observed
deletions, red circles denote imputed deletions, and red squared denote imputed non-deletions.  
A solid vertical line separates the marker blocks.

A separate plot is produced for each chromosome. 

The marker order is assumed to be the same as the column order of \code{x}.


}
\author{
Andrew George
}
\seealso{
\code{\link{idmarkerblocks}}
}
\examples{
# generate map for two chromosomes that have 50 equally distance marker loci
map <- sim.dmmap(rep(100,2), 50)

# generate deletion data
Slist <-  sim.dmdata( numlines = 250, plines = 0.2, map = map, 
    Enumdel = 4, p.missing = 0.02,  seed = 1)


# produce deletion map plot for each chromosome when data have no
# missing marker genotypes
dmplot(Slist[["nomissing"]], main="No missing", keep.chrm=1)

# produce deletion map plot for each chromosome when data have 
# missing marker genotypes. Here, no blocking structure is given.
dmplot(Slist[["missing"]], main="Missing", keep.chrm=1)

}

