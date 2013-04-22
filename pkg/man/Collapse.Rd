\name{Collapse}
\alias{Collapse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Collapse Redundant Marker Loci
}
\description{
Redundant marker loci are collapsed to a single marker locus. 
}
\usage{
Collapse(dmat= NULL)
}
\arguments{
  \item{dmat}{
An object of class \code{\link{dmdata}}
}
}
\details{
Equivalent marker loci are collapsed to a single unique marker locus. Marker loci that contain missing 
genotypes will not be collapsed.
}
\value{
A list is returned, containing two elements; data and collapsed. The \code{data} element contains the collapsed data 
and is of type \code{\link{dmdata}}. The \code{collapsed} element is a list of length the number of groups 
of marker loci that have been collapsed into a single marker locus. Each element in the list is indexed by the 
unique marker locus upon which the redundant marker loci have been collapsed.  
}
\author{
Andrew W George and James Totterdell
}
\seealso{
\code{\link{CleanData}}
}
\examples{
# generate map for a single chromosome with 20 equally distant marker loci
map <- sim.dmmap(100,20)



# generate deletion data
Slist <-  sim.dmdata( numlines = 50, plines = 0.2, map = map, 
    Enumdel = 4, p.missing = 0.02,  seed = 1)


# remove redundant marker loci when there is no missing marker data
a <- Collapse(Slist[["nomissing"]])

# a$data contains 50 rows and 11 columns
dim(a$data)

# now removed redundant marker loci where there are missing marker data
a <- Collapse(Slist[["missing"]])

# a$data now contains 50 rows and 19 columns because marker loci that contain NAs cannot be collapsed. 
dim(a$data)

}
