\name{impute}
\alias{impute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Impute Missing Marker Genotypes
}
\description{
This function imputes any missing marker genotypes.  
}
\usage{
impute(x, uniform=TRUE)
}
\arguments{
  \item{x}{
An object of class \code{\link{dmdata}} or \code{\link{dmdatachrm}} which contains 
the genotype data in 0/1 form and any 
missing marker genotypes.
}
    \item{uniform}{
A logical variable that when true, specifies that the imputation of any missing marker genotypes 
will be uniform, otherwise, see below.
}
}
\details{
\code{impute} is used in \code{\link{dmapping}} routine for imputing missing genotypes. 

When \code{uniform} is \code{TRUE}, the imputation of the missing genotypes (0,1) is uniform. 

When \code{uniform} is \code{FALSE}, the imputation of a missing genotype is dependent upon its neighbours. 
If the flanking genotypes are both 0, then the missing genotype is assigned 0 with probability 0.95. 
If one of the flanking genotypes is 1, then the missing genotype is assigned 0 with probability 0.5.
If the flanking genotypes are both 1, then the missing genotype is assigned 0 with probability 0.05. 
}
\value{
An object of class \code{\link{dmdata}} or \code{\link{dmdatachrm}} is returned that has the 
missing genotypes imputed. 
}
\author{
Andrew W. George
}
\seealso{
\code{\link{dmapping}}
}
\examples{
# generate map with 10 chromosomes, each of length 
# 100 cM and each with 50  equally spaced marker loci. 
map <- sim.dmmap(rep(100,10), 50)

# generate deletion data
Slist <-  sim.dmdata( numlines = 50, plines = 0.2, map = map, 
    Enumdel = 4, p.missing = 0.02,  seed = 1)

# class of object that does not contain any missing data
class(Slist[["missing"]])


# impute missing genotypes where 0,1 are equally likely
impute(Slist[["missing"]])

# impute missing genotypes conditional on flanking genotypes. 
impute(Slist[["missing"]], uniform=FALSE)

}

