\name{SimDeletions}
\alias{SimDeletions}
\alias{simulate}
\alias{simulate deletion data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation
}
\description{
'SimDeletions' returns a matrix of genotype data where the rows correspond to 
plants and the columns correspond to tightly linked marker loci. A 
column (i.e., marker locus) can contain genotypes 0/1/NA where 0 is no 
deletion, 1 is a deletion, and NA is missing.   
}
\usage{
SimDeletions(numlines = NULL, plines = 0.5, nummarkers = NULL, labels = NULL, Enumdel = 4, p.missing = NULL, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numlines}{
number of inbred lines/plants.
}
  \item{plines}{
the probability of a plant carrying a deletion. 
}
  \item{nummarkers}{
the number of marker loci spanning the genomic region of interest.
}
  \item{labels}{
a character vector containing the names of the marker loci. 
}
  \item{Enumdel}{
the average length (in terms of number of marker loci) for a deletion.
}
  \item{p.missing}{
the probability of a genotype being unobserved.
}
  \item{seed}{
an integer number that if set, will produce reproducable replicates. 
}
}
\details{
Data are generated assuming the marker loci are in the true order 
and that a plant can only contain a single deletion (but that deletion may 
span several marker loci).  The length of a deletion is random and follows a 
poisson distribution with mean Enumdel. If p.missing is set, then the 
data will contain missing genotypes (coded NA). 
}
\value{
a list is returned that contains only a single element if p.missing is NULL 
or two elements if p.missing is set. The elements are:
\item{nomissing}{is a matrix of genotype data 0/1 where 0 is no deletion and
1 is a deletion. The number of rows corresponds to the number of plants 
in the simulated population and the number of columns is the number of 
genetic marker loci. This will be the only element in output list if 
p.missing is NULL.}
\item{missing}{is a matrix of genotype data 0/1 that  
is equivalent to "nomissing" but where marker genotypes have 
been replaced, randomly, and set to NA.} 
}
\author{
Andrew W. George
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{plot}}
}
\examples{

# To generate a inbred population of 500 plants 
# where genotype data are collected on 20 marker loci,
# the probability of a plant having a genomic deletion is 0.2, 
# the average length of a deletion is 8 marker loci, and 
# the probability of a marker genotype being missing is 0.02,
# use:
S <- SimDeletions( numlines = 500, plines = 0.2, nummarkers = 20,  
    Enumdel = 8, p.missing = 0.02,  seed = NULL)

# print contents of matrix
print(S$nomissing)

# print contents of matrix
print(S$missing)
}
 
