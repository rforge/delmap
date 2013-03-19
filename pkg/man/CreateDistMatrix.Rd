\name{CreateDistMatrix}
\alias{CreateDistMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distance Measure Computation
}
\description{
This function computes the Manhattan distance matrix from 0/1 deletion data. It is assumed the
deletion data does not contain any missing genotypes. 
}
\usage{
CreateDistMatrix(dmat = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmat}{
a \code{delmap.data} object or matrix that 
contains only 0's and 1's.
}
}
\details{
Calculates the Manhattan distance between all pairs of marker loci (columns). The function 
will produce an error if the data contains missing values. 
}
\value{
Returns an object of class 'dist'.
}
\author{
Andrew W. George
}

\examples{
# deletion data where 1 is a deletion
d <- matrix(data=c(0,0,1,1, 
                        0,0,0,1, 
                        1,0,0,0, 
                        0,1,0,0,  
                        1,1,1,0), nrow=5)
# convert  into deldata.map  object
dmap <- as.delmap(d)

# calculate distance matrix for deletion data
CreateDistMatrix(dmat=dmap)

}

