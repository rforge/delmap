\name{ImputeMissingGeno}
\alias{ImputeMissingGeno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Impute Missing Marker Genotypes
}
\description{
This function imputes, uniformly,  the missing marker genotypes.  
}
\usage{
ImputeMissingGeno(dmat = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmat}{
An object of class \code{\link{delmap.data}} which contains the genotype data in 0/1 form and contains missing genotypes.
}
}
\details{
Missing genotypes are imputed, assuming 0,1 occurs with equal frequency. 
}
\value{
An object of class \code{\link{delmap.data}} is returned that has the missing genotyeps imputed. 
}
\author{
Andrew W. George
}
\examples{
# example to be done
# generate deletion data
Slist <-  SimDeletions( numlines = 50, plines = 0.2, nummarkers = 20,
    Enumdel = 4, p.missing = 0.02,  seed = 1)

# class of object that does not contain any missing data
class(Slist[["missing"]])


# impute missing genotypes
ImputeMissingGeno(Slist[["missing"]])



}

