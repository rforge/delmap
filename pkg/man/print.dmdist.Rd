\name{print.dmdist}
\alias{print.dmdist}
\title{
Print Deletion Data
}
\description{
'print' prints the deletion data.
}
\usage{
\method{print}{dmdist}(x, ...)
}
\arguments{
  \item{x}{
an object of class \code{\link{dmdist}}.
}
  \item{...}{
further arguments passed to or from other methods.
}
}
\details{
\code{print.dmdist} is a customized print function for printing the contents of the \code{dmdist} object. 
The contents are printed in tabular form. The first column is the number of marker blocks present in the reference 
block structure. The second column is the number of blocks present in the observed block structure.   
The third column is the minimum rearrangement distance summed across the observed blocks. The fourth 
column is the total number of wrongly placed marker loci. 
}
\author{
Andrew W. George
}
\seealso{
\code{\link{RearrangeDist}}
}
\examples{
# need an example !!!!!!!!!
}

