\name{delmap-package}
\alias{delmap-package}
\alias{delmap}
\docType{package}
\title{
Deletion Mapping
}
\description{
A package for using deletion data to fine map genetic marker loci. It 
is assumed the marker loci are (almost) completely linked and unable to be 
ordered through linkage. It is also assumed that a plant will contain, at 
most, a single deletion that spans multiple marker loci. The deletions 
are allowed to vary in size and there may be unobserved marker data. The 
package contains functions to find the best ordering of the marker loci, 
to visualize the deletion map, and to impute missing marker genotypes.   
}
\details{
\tabular{ll}{
Package: \tab delmap\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-11-12\cr
License: \tab GPL (>=2)\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Author:  Andrew W. George
Maintainer: Andrew W. George <andrew.george@csiro.au>
}
\keyword{ package }
\examples{
S <- SimDeletions(numlines = 500, plines = 0.2, nummarkers = 20,
    Enumdel = 8, p.missing = 0.02,  seed = NULL)
}
