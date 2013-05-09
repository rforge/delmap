\name{RearrangeDist}
\alias{RearrangeDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rearrangement Distance
}
\description{
A multivariate metric to measure the "distance" between two deletion mapping orders. 
}
\usage{
RearrangeDist(refblockstr=NULL, refmrkord=NULL, obsblockstr=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{refblockstr}{
a integer vector for the reference block structure. It is a named (with the marker labels) vector. 
}
  \item{refmrkord}{
a list object of length the number of unique blocks in the reference ordering. Each element is an integer vector 
giving the marker ordering with the block. It is a named (with the marker labels) vector. 
}
  \item{obsblockstr}{
a integer vector for the observed/realized block structure. It is a named (with the marker labels) vector. 
}
}
\details{

\code{RearrangeDist} is a multivariate metric for measuring the distance between a reference order
and an observed order. It is based on the L1 rearrangement distance of Amir et al. (2009) Journal of Computer and Systems 
Sciences. Amir's rerrangement distance is a  measure for measuring the distance between two strings. However, 
comparing two deletion mappign orders is not as simple as comparing two strings. 
Two orders may differ in their block size, number of blocks, which markers are assigned to which block, and how the markers 
are ordered within a block.  Hence, \code{RearrangeDist} reports a number of measures of distance.   

To describe \code{RearrangeDist},  consider the following example. 

Suppose the reference block structure is
\tabular{cccccccccc}{
M1 \tab M2 \tab M3 \tab M4 \tab M5 \tab \tab M6 \tab M7 \tab M8 \tab M9 \cr
1  \tab 1  \tab  1 \tab 1  \tab 1  \tab \tab 2  \tab 2  \tab 2  \tab 2  \cr
}
Here, markers M1 to  M5 have overlapping deletions and hence belong to the same block structure. 
Similarly, markers M6 to M9 have overlapping deletions (but not with M1 to M5) and belong to the same block 
structure. 

Let the marker ordering for the reference order be 
\tabular{cccccccccc}{
M1 \tab M2 \tab M3 \tab M4 \tab M5 \tab \tab M6 \tab M7 \tab M8 \tab M9 \cr
1  \tab 1  \tab  2 \tab 2  \tab 3  \tab \tab 1  \tab 1  \tab 2  \tab 2  \cr}
Here, the deletion data for M1 and M2 are identical, hence they are both indexed by 1. Markers 
M3 and M4 also have identical data in this example, and are both indexed by 2. Marker M5 is different to the other 
markers in the block so is assigned index 3. Similarly, for the markers in the second block. Markers M6 and M7 are identical and  
are given index 1. Markers M8 and M9 are identical and given index 2.  

Also, let the block structure for the observed marker order be 
\tabular{ccccccccccc}{
M1 \tab M6 \tab M5 \tab M4 \tab \tab M2 \tab M7 \tab \tab M3 \tab M9 \tab M8  \cr
1  \tab 1  \tab 1  \tab 1   \tab \tab 2 \tab 2  \tab \tab 3  \tab 3  \tab 3   \cr
}
By running \code{RearrangeDist} on these data, it would output 
\itemize{
\item \code{nrefblocks} the number of blocks in the reference order which is 2
\item \code{nobsblocks} the number of blocks in the observed order which is  3
\item \code{mindist} the total minimum rearrangement distance which is 2+0+0=2 
\item \code{minwrong} the total minimum number of wrongly placed markers is 1+1+1=3
}
See the vignette for details on how \code{mindist} is calculated. 
}
\value{
A list with elements \code{nrefblocks} for the number of blocks in the reference ordering, 
\code{nobsblocks} for the number of blocks in the obsered marker ordering, 
\code{mindist} for the minimum total rearrangement distance, and
\code{minwrong} for the minumum total number of misplaced marker loci. A misplaced marker locus is a marker locus 
belong to the wrong block of marker loci. 
}
\author{
Andrew W George
}
\note{
In calculating our distance measure, we do not require a list object giving the marker ordering within 
each block for the observed order. Only the blocking structure for the observed order is required.  
}
\seealso{
\code{\link{idmarkerblocks}}, \code{\link{idmarkersuborder}}, \code{\link{dmapping}}
}
\examples{
# reference block structure
bl <- c(rep(1, 5), rep(2, 4))
names(bl) <-  paste("M", 1:9, sep="")

# reference marker ordering
mr <- list()
mr[[1]] <- c(1,1,2,2,3)
mr[[2]] <- c(1,1,2,2)
names(mr[[1]]) <- paste("M", 1:5, sep="")
names(mr[[2]]) <- paste("M", 6:9, sep="")

# observed block structure for which we want to calculate a distance
ob <- c(1,1,1,1,2,2,3,3,3)
names(ob) <- c("M1", "M6", "M5", "M4", "M2", "M7", "M3", "M9", "M8")

# calculate rearrangement distance between reference and observed order 
RearrangeDist(refblockstr=bl, refmrkord=mr, obsblockstr=ob)


}
