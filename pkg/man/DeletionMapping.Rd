\name{DeletionMapping}
\alias{DeletionMapping}
\title{ Deletion Mapping
}
\description{
\code{DeletionMapping} orders the marker loci through the pattern of overlapping observed genomic deletions. 
}
\usage{
DeletionMapping(dmap=NULL, niterates=100, nwithin=100, temp=1000, psampled=0.1, method="concorde",...)
}
\arguments{
  \item{dmap}{
An object of class \code{\link{delmap.data}}
}
  \item{niterates}{
specifies the number of iterations of the simulated annealing algrorithm to perform
}
  \item{nwithin}{
specifies the number of iterations of the simulated annealing algorithm to perform 
for a given temperature
}
 \item{temp}{
starting temperature for the simulated annealing algorithm.
}
 \item{psampled}{
the proportion of missing genotypes to impute with new values at each iteration of the simulated 
annealing algorithm. 
}
 \item{method}{method to solve the TSP (default: nearest insertion
          algorithm; see details).
}
 \item{...}{other parameters that are to be passed to the TSP solving function \code{solve_TSP}.
}
}
\details{
Simulated annealing ...
}
\value{
An object of class \code{\link{delmap.data}}.
}
\author{
Andrew W George
}
\seealso{
\code{\link{Collapse}}
}
\examples{
d <-   SimDeletions( numlines = 50, plines = 0.2, nummarkers = 50,
    Enumdel = 8, p.missing = 0.1,  seed = 1)

dc <- CleanData(d[["missing"]])

permdc <- PermuteCols(dc)
a <- DeletionMapping(dmap=permdc, niterates=50, nwithin=50, temp=50, psampled=0.2, method="2-opt")


# to run the exact TSP solver concorde, first specify the absolute path of where it 
# has been installed
# concorde_path(path="/home/cmis/geo047/Condore")
#a2 <- DeletionMapping(dmap=permdc, niterates=50, nwithin=50, temp=1000, psampled=0.2, method="concorde")
}
