\name{ImputeMissingGeno}
\alias{ImputeMissingGeno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Impute Missing Marker Genotypes
}
\description{
This function imputes any missing marker genotypes.  
}
\usage{
ImputeMissingGeno(dmat = NULL, uniform=TRUE)
}
\arguments{
  \item{dmat}{
An object of class \code{\link{delmap.data}} which contains the genotype data in 0/1 form and any 
missing marker genotypes.
}
    \item{uniform}{
A logical variable that when true, specifies that the imputation of any missing marker genotypes 
will be uniform, otherwise, see below.
}
}
\details{
\code{ImputeMissingGeno} is used in \code{\link{DeletionMapping}} routine for imputing missing genotypes. 

When \code{uniform} is \code{TRUE}, the imputation of the missing genotypes (0,1) is uniform. 

When \code{uniform} is \code{FALSE}, the imputation of a missing genotype is dependent upon its neighbours. 
If the flanking genotypes are both 0, then the missing genotype is assigned 0 with probability 0.95. 
If one of the flanking genotypes is 1, then the missing genotype is assigned 0 with probability 0.5.
If the flanking genotypes are both 1, then the missing genotype is assigned 0 with probability 0.05. 
}
\value{
An object of class \code{\link{delmap.data}} is returned that has the missing genotyeps imputed. 
}
\author{
Andrew W. George
}
\seealso{
\code{\link{DeletionMapping}}
}
\examples{
# example to be done
# generate deletion data
Slist <-  SimDeletions( numlines = 50, plines = 0.2, nummarkers = 20,
    Enumdel = 4, p.missing = 0.02,  seed = 1)

# class of object that does not contain any missing data
class(Slist[["missing"]])


# impute missing genotypes where 0,1 are equally likely
ImputeMissingGeno(Slist[["missing"]])[,]

# impute missing genotypes conditional on flanking genotypes. 
ImputeMissingGeno(Slist[["missing"]], uniform=FALSE)[,]

}

