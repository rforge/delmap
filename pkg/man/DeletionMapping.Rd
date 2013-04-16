\name{DeletionMapping}
\alias{DeletionMapping}
\title{ Deletion Mapping
}
\description{
\code{DeletionMapping} orders the marker loci through the pattern of overlapping observed genomic deletions. 
}
\usage{
DeletionMapping(dmap=NULL, niterates=100, nwithin=100, cooling=0.99, psampled=0.1, method="concorde",
                refblockstr=NULL, refmrkord=NULL, ...)
}
\arguments{
  \item{dmap}{
An object of class \code{\link{dmdata}} that contains coded 0, 1 marker genotypes and missing 
marker genotypes (NAs). 
}
  \item{niterates}{
specifies the number of iterations of the simulated annealing algrorithm to perform
}
  \item{nwithin}{
specifies the number of iterations of the simulated annealing algorithm to perform 
for a given temperature
}
 \item{cooling}{
coefficient (between 0 and 1) for geometric cooling of temperature for the simulated annealing algorithm.
}
 \item{psampled}{
the proportion of missing genotypes to impute with new values at each iteration of the simulated 
annealing algorithm. 
}
 \item{method}{method to solve the TSP (default: nearest insertion
          algorithm; see details).
}
 \item{refblockstr}{integer vector with true block assignments if known.
}
 \item{refmrkord}{list of length the number of blocks. Each element is an integer vector giving the order of the markers within the block.}  
 \item{...}{other parameters that are to be passed to the TSP solving function \code{solve_TSP}.
}
}
\details{
Simulated annealing ...
}
\value{
An object of class \code{\link{dmdata}}.
}
\author{
Andrew W George
}
\seealso{
\code{\link{Collapse}}
}
\examples{

# simulate deletion data with 50 marker loci collected on 50 inbred lines
dfull <-   SimDeletions( numlines = 50, plines = 0.2, nummarkers = 50,
    Enumdel = 8, p.missing = 0.1,  seed = 1)

# clean deletion data by removing noninformative rows/columns
dclean <- CleanData(dfull[["missing"]])

# obtain blocking structure  for the simulated true data
mrks <- colnames(dclean)
cindx <- match(mrks, colnames(dfull[["nomissing"]]))
dnoNA <- as.dmdata(dfull[["nomissing"]][,cindx])
refblockstr <- IdentifyMarkerBlocks(dnoNA)  ## simulated true block structure
refmrkord   <- IdentifyMarkerOrd(dnoNA)$orders  ## simulated true marker structure


# permute dclean so simulated data not in true order
dperm <- PermuteCols(dclean)
a <- DeletionMapping(dmap=dperm, niterates=50, nwithin=5, cooling=0.95, psampled=0.2, method="2-opt", 
                     refblockstr=refblockstr, refmrkord=refmrkord)

# colour plot of deletion map
plot(a$bestmap)




# to run the exact TSP solver concorde, first specify the absolute path of where it 
# has been installed
# concorde_path(path="/home/cmis/geo047/Condore")
#a2 <- DeletionMapping(dmap=permdc, niterates=50, nwithin=50, temp=1000, psampled=0.2, method="concorde")
}
