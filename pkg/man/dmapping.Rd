\name{dmapping}
\alias{dmapping}
\title{ Deletion Mapping
}
\description{
\code{dmapping} orders the marker loci through the pattern of overlapping observed genomic deletions. 
}
\usage{
dmapping(x=NULL, chrm=NULL, niterates=100, nwithin=100, cooling=0.99, psampled=0.1, method="concorde",
                refblockstr=NULL, refmrkord=NULL, ...)
}
\arguments{
  \item{x}{
An object of class \code{\link{dmdata}} or \code{\link{dmdatachrm}} that contains the deletion data.
}
  \item{chrm}{
an integer vector indexing which chromosomes are to be analyzed.
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
A list object with elements \code{bestx}, \code{bestt},
\code{besttvec},  and if reference ordering and block structure is provided,
\code{bestdist}.

Here, \itemize{
\item \code{bestx} is of object type \code{\link{dmdatachrm}}  and contains the best
deletion map with imputed genotypes.
\item \code{bestt} is a \code{numeric} value containing the (best) tour length for
\code{bestx}.
\item \code{besttvec} is a \code{numeric} vector containing the best tour lengths at each
iteration of the simulated annealing algorithm.
\item \code{bestdist} is a \code{\link{dmdist}} object for the (best) rearrangement distance
for \code{bestx}.
}

}
\author{
Andrew W George
}
\seealso{
\code{\link{dmappingchrm}}
}
\examples{
# generate map for a single chromosome with 20 equally distant marker loci
map <- sim.dmmap(100,20)

# simulate deletion data with 50 marker loci collected on 50 inbred lines
dfull <-   sim.dmdata( numlines = 50, plines = 0.2, map = map, 
    Enumdel = 8, p.missing = 0.1,  seed = 1)

# clean deletion data by removing noninformative rows/columns
dclean <- cleandata(dfull[["missing"]])

# obtain blocking structure  for the simulated true data
#mrks <- colnames(dclean)
#cindx <- match(mrks, colnames(dfull[["nomissing"]]))
#dnoNA <- as.dmdata(dfull[["nomissing"]][,cindx])
#refblockstr <- IdentifyMarkerBlocks(dnoNA)  ## simulated true block structure
#refmrkord   <- IdentifyMarkerOrd(dnoNA)$orders  ## simulated true marker structure

# permute dclean so simulated data not in true order
dperm <- permutecols(dclean)
a <- dmapping(x=dperm, niterates=50, nwithin=50, cooling=0.95, psampled=0.2, method="2-opt", 
                     refblockstr=NULL, refmrkord=NULL)

#x=dperm[[1]]; niterates=50; nwithin=5; cooling=0.95; psampled=0.2; method="2-opt"


# colour plot of deletion map
dmplot(a$bestx, cx=0.75)




# to run the exact TSP solver concorde, first specify the absolute path of where it 
# has been installed
# concorde_path(path="/home/cmis/geo047/Condore")
#a2 <- dmapping(dmap=permdc, niterates=50, nwithin=50, temp=1000, psampled=0.2, method="concorde")
}
