\name{CreateDistMatrix}
\alias{CreateDistMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distance Measure Computation
}
\description{
This function computes a distance matrix from 0/1 deletion data. It is assumed the
deletion data does not contain any missing genotypes. 
}
\usage{
CreateDistMatrix(mat = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
a matrix or an object that can be converted to be a matrix that 
contains only 0's and 1's.
}
}
\details{
'CreateDistMatrix' is a function for calculating a distance matrix between all pairs of loci 
from deletion data. This matrix is foundational to being able to order marker loci based on deletions. 
Our underlying assumption is that loci that are very close to each other will 
share a greater proportion of deletions across lines than those loci that are far apart. We
expect this assumption to be 
valid if the genomic window being spanned by the marker loci is very small.

Our distance measure is based on the number of deletions in common between two loci. 
Suppose the gentoype data  are
\tabular{cccc}{
M1 M2 M3 M4 \cr
0    0   1   1  \cr
0    0   0   1  \cr
1    0   0   0  \cr
0    1   0   0  \cr
1    1   1   0  
}
where we have data on five inbred plants on four marker loci, and 1 is a deletion.

First, the number of lines that share a deletion between two loci is calculated, giving
\tabular{cccc}{
    2  1  1  0 \cr
    1  2  1  0 \cr
    1  1  2  1 \cr
    0  0  1  2 
}

Second, we subtract this matrix from its maximum value giving
\tabular{cccc}{
  0    1    1    2 \cr
  1    0    1    2 \cr
  1    1    0    1 \cr
  2    2    1    0
}

Third, the diagonal elements are set to a large value (say 100) giving the distance matrix
\tabular{rrrr}{
  100    1    1    2 \cr
  1    100    1    2 \cr
  1    1    100    1 \cr
  2    2    1    100 
}

We can then use this distance matrix in the Travelling Salesman Problem to find the "best" ordering of the marker loci. 

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns an object of class 'dist'.
}
\author{
Andrew W. George
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index) for the standard data sets.

# deletion data where 1 is a deletion
dd <- matrix(data=c(0,0,1,1, 
                        0,0,0,1, 
                        1,0,0,0, 
                        0,1,0,0,  
                        1,1,1,0), nrow=5)
# calculate distance matrix for deletion data
CreateDistMatrix(mat=dd)

}

