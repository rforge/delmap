\name{CleanData}
\alias{CleanData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clean Data in Readiness for Analysis
}
\description{
\code{CleanData} removes those rows and columns in the deletion data that do not carry any deletions.
}
\usage{
CleanData(dmat= NULL, ignoreNA=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmat}{
An object of class \code{\link{delmap.data}}
}
  \item{ignoreNA}{
if set to TRUE, NA's are ignored when determining if a row/column should be removed. If set to FALSE, 
rows/columns that carry more than one NA will be retained. 
}
}
\details{
Rows and columns of genotype data that have  no deletions are removed.  By default, NA's are ignored. If 
\code{ignoreNA} is set to FALSE, then rows/columns with NA's will 
not be removed. 

However, if a row or column has only a single NA, then it is still removed, even if 
\code{ignoreNA} is set to FALSE. 

}
\value{
An object of class \code{\link{delmap.data}}.
}
\author{
Andrew W George
}
\seealso{
\code{\link{ReadData}}, \code{\link{ImputeMissingGeno}}, \code{\link{Collapse}}
}
\examples{
# generate deletion data
Slist <-  SimDeletions( numlines = 50, plines = 0.2, nummarkers = 20,
    Enumdel = 4, p.missing = 0.02,  seed = 1)

# Number of marker loci before beign cleaned
nrow(Slist[["missing"]])

# Clean data but where missing genotypes are ignored.
a <- CleanData(Slist[["missing"]])
a[,]  # print data matrix
nrow(a)

# Clean data but where colums/rows with missing genotypes are retained
# for later use. 
b <- CleanData(Slist[["missing"]],ignoreNA=FALSE)
b[,] #print data matrix
nrow(b)

}
