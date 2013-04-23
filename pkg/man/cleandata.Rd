\name{cleandata}
\alias{cleandata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clean Data 
}
\description{
\code{cleandata} removes those rows and columns that do not carry any deletions. See \code{details} for 
impact of NAs.
}
\usage{
cleandata(dmdat= NULL, ignoreNA=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmdat}{
An object of class \code{\link{dmdata}}
}
  \item{ignoreNA}{
A logical variable. When TRUE, NA's are ignored and a row/column is only retained if it carrys at least 
one deletion. When FALSE, a column will be removed only if it has no deletions or NAs. A row will be 
removed if it has no deletions. Rows are removed first, followed by columns.  
}
}
\details{
When \code{ignoreNA} is TRUE, row and columns that do not carry any deletions are removed. 

When \code{ignoreNA} is FALSE, first rows that do not carry any deletions are removed. Then, of the remaining
deletion data, those columns that carry no deletions and have no NAs are also removed.  
Here, the assumption is that if a plant(row) has no recorded 
deletions, then it is unlikely that the NA's are deletions.  
}
\value{
An object of class \code{\link{dmdata}}.
}
\author{
Andrew W George
}
\seealso{
\code{\link{read.dmdata}}, \code{\link{ImputeMissingGeno}}, \code{\link{Collapse}}
}
\examples{
# generate map for a single chromosome with 20 equally distant marker loci
map <- sim.dmmap(100,50)

# generate deletion data
Slist <-  sim.dmdata( numlines = 50, plines = 0.2, map = map,
    Enumdel = 4, p.missing = 0.05,  seed = 1)

# Number of marker loci before being cleaned
nrow(Slist[["missing"]])


# Clean data but where NAs affect which rows/columns are removed
a <- cleandata(Slist[["missing"]])
a[,] #print data matrix
nrow(a)

}