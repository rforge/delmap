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
if set to TRUE, missing genotypes will be ignored. If set to FALSE, rows/columns with missing genotypes will be retained.
}
}
\details{
Rows and columns of genotype data that have  no deletions can be removed from the analysis. This has the advantage of 
reducing the dimension of the problem. 
}
\value{
An object of class \code{\link{delmap.data}}.
}
\author{
Andrew W George
}
\note{
By default, missing genoyptes are ignored, which will lead to rows/columns being deleted if they do not carry 
any deletions. However, it may be more approriate to set \code{ignoreNA=FALSE}. This will ensure rows/columns that 
have missing genotypes will be retained. Imputation can then be used later to impute the missing genotypes. 
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
CleanData(Slist[["missing"]])
nrow(CleanData(Slist[["missing"]]))

# Clean data but where colums/rows with missing genotypes are retained
# for later use. 
nrow(CleanData(Slist[["missing"]],ignoreNA=FALSE))
}
